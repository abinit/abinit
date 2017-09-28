!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_effpot_mpi
!!
!! NAME
!! m_effpot_mpi
!!
!! FUNCTION
!! Module for using the paralelisation of effective potential
!! Container type is defined, and destruction, print subroutines
!! This module is still experimental
!!
!! COPYRIGHT
!! Copyright (C) 2010-2017 ABINIT group (AM)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_effpot_mpi

 use defs_basis
 use m_xmpi
 use m_errors
 use m_profiling_abi

 implicit none

 public :: effpot_mpi_init
 public :: effpot_mpi_free
!!***

!!****t* defs_abitypes/effpot_mpi_type
!! NAME
!! effective_potential_type
!!
!! FUNCTION
!! datatype to set the parallelisation
!!
!! SOURCE

 type, public :: effpot_mpi_type

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This is for the parallelisation over the supercell
   integer :: comm
   ! local Communicator over all processors treating the same cell

   integer :: my_rank
   ! local Communicator over all processors treating the same cell

   integer my_ncell
   ! Number of cell treat by current proc
   
   integer my_nrpt
   ! Number of rpt/coefficient treat by current proc

   integer,allocatable :: my_cells(:)
   ! my_cells(my_ncell)
   ! Number of each cells in the supercell treat by this CPU

   integer,allocatable :: my_index_cells(:,:)
   ! my_cells(my_ncell,3)
   ! indexes of the cells in the supercell treat by this CPU

   integer,allocatable :: my_rpt(:)
   ! my_rpt(my_nrpt)
   ! List of rpt treat by this CPU

 end type effpot_mpi_type
!!***

CONTAINS  !===========================================================================================

!****f* m_effective_potential/effpot_mpi_init
!!
!! NAME
!! effpot_mpi_init
!!
!! FUNCTION
!! deallocate all dynamic memory for mpi of supercell
!!
!! INPUTS
!! cell_number(3) = size of the supercell, 3 3 3 for example
!! ndiv = number of division to consider. For example if ndiv==2,
!!        the mpi will be set over the 2 lvl of parallelisation,
!!        cell and nrpt by considering nrpt / 2 for each CPU 
!!        (still experimental, only ndiv=1 available)
!! nrpt = number of cell to be parallelised for the IFC
!! comm = MPI communicator
!!
!! OUTPUT
!! effpot_mpi<type(effpot_mpi_type)()> = effpot_mpi datatype
!! PARENTS
!!      m_effective_potential
!!
!! CHILDREN
!!
!! SOURCE

subroutine effpot_mpi_init(cell_number,effpot_mpi,ndiv,nrpt,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'effpot_mpi_init'
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm,ndiv,nrpt
!array
 type(effpot_mpi_type), intent(inout) :: effpot_mpi
 integer,intent(in) :: cell_number(3)
!Local variables-------------------------------
!scalars
 integer :: i1,i2,i3,icell,ii,irpt
 integer :: my_rank,ncell_alone,ncell,nproc
 integer :: npcell,nprpt,virt_rank
 integer :: master = zero
 logical :: iam_master = .false.
 character(len=500) :: msg
!array
 integer,allocatable :: rpt_list(:)

! ***********************************************************************

!Set the number of cell in the supercell
 ncell = product(cell_number(:))


 if (any(cell_number <= 0).or.ncell<=0) then
   write(msg,'(a,a)')' No supercell found for setting'
   MSG_ERROR(msg)
 end if

!MPI variables
 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 iam_master = (my_rank == master)


!Set the number of cpu for each level
 npcell = ncell / nproc * ndiv
 nprpt  = nrpt  / ndiv

!Do some checks
 if(mod(nrpt,npcell) /= 0.or.mod(ncell,npcell) /=0)then
!   write(msg,'(2a,2I0)')' Chose another number of CPU ',ncell,nrpt
!   MSG_ERROR(msg)   
 end if

 call effpot_mpi_free(effpot_mpi)

 effpot_mpi%comm = comm
 effpot_mpi%my_rank = my_rank

!Determine the number of cell for each CPU
 ncell_alone = mod(ncell,nproc)

!TREAT CELL
 effpot_mpi%my_ncell = aint(real(ncell,sp)/(nproc))

 if(my_rank >= (nproc-ncell_alone)) then
   effpot_mpi%my_ncell = effpot_mpi%my_ncell  + 1
 end if

 if(ndiv>1) effpot_mpi%my_ncell = npcell

!Allocation of array
 ABI_ALLOCATE(effpot_mpi%my_cells,(effpot_mpi%my_ncell))
 ABI_ALLOCATE(effpot_mpi%my_index_cells,(effpot_mpi%my_ncell,3))
 effpot_mpi%my_cells = zero
 effpot_mpi%my_index_cells = zero

 virt_rank = aint(real(my_rank,sp)/(ndiv))
 
 do icell=1,effpot_mpi%my_ncell
   if(virt_rank >= (nproc-ncell_alone))then
     effpot_mpi%my_cells(icell)=(aint(real(ncell,sp)/nproc))*(virt_rank)+&
&                              (virt_rank - (nproc-ncell_alone)) + icell
   else
     effpot_mpi%my_cells(icell)=(effpot_mpi%my_ncell)*(virt_rank)  + icell
  end if
 end do

 icell = 0
 ii = 0
 do i1 = 1,cell_number(1)
   do i2 = 1,cell_number(2)
     do i3 = 1,cell_number(3)
       ii = ii +1
       if(any(effpot_mpi%my_cells==ii))then
         icell=icell+1
         effpot_mpi%my_index_cells(icell,1) = i1;
         effpot_mpi%my_index_cells(icell,2) = i2;
         effpot_mpi%my_index_cells(icell,3) = i3;
       end if
     end do
   end do
 end do

!TREAT RPT Not yet parallelised 
 ABI_ALLOCATE(rpt_list,(nproc))
 
 irpt = 1

 do while (irpt < nproc)
   do ii=1,ndiv       
     rpt_list(irpt+ii-1) = ii-1
   end do
   irpt = irpt + ndiv
 end do

 virt_rank = rpt_list(my_rank+1)
 
 effpot_mpi%my_nrpt = nprpt
 ABI_ALLOCATE(effpot_mpi%my_rpt,(effpot_mpi%my_nrpt))
 do irpt=1,effpot_mpi%my_nrpt
!AM_EXPERIMENTAL
!   if(virt_rank >= (nproc-ncell_alone))then
!    effpot_mpi%my_rpt(irpt)=(aint(real(nrpt,sp)/nproc))*(virt_rank)+&
!&                              (virt_rank - (nproc-ncell_alone)) + irpt
!   else
!     effpot_mpi%my_rpt(irpt)=(effpot_mpi%my_nrpt)*(virt_rank)  + irpt
!   end if
     effpot_mpi%my_rpt(irpt)= irpt
!AM_EXPERIMENTAL
 end do

 ABI_DEALLOCATE(rpt_list)

end subroutine effpot_mpi_init
!!***

!****f* m_effective_potential/effpot_mpi_free
!!
!! NAME
!! effpot_mpi_free
!!
!! FUNCTION
!! deallocate all dynamic memory for mpi
!!
!! INPUTS
!! effpot_mpi<type(effpot_mpi_type)()> = effpot_mpi datatype
!!
!! OUTPUT
!!
!! PARENTS
!!      m_effective_potential,m_effpot_mpi
!!
!! CHILDREN
!!
!! SOURCE

subroutine effpot_mpi_free(effpot_mpi)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'effpot_mpi_free'
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
!array
 type(effpot_mpi_type), intent(inout) :: effpot_mpi

!Local variables-------------------------------
!scalars
!array

! *************************************************************************

 effpot_mpi%my_ncell = 0
 effpot_mpi%my_nrpt = 0

 if (allocated(effpot_mpi%my_cells)) then
   effpot_mpi%my_cells(:) = zero
   ABI_DEALLOCATE(effpot_mpi%my_cells)
 end if

 if (allocated(effpot_mpi%my_index_cells)) then
   effpot_mpi%my_index_cells(:,:) = zero
   ABI_DEALLOCATE(effpot_mpi%my_index_cells)
 end if

 if (allocated(effpot_mpi%my_rpt)) then
   effpot_mpi%my_rpt(:) = zero
   ABI_DEALLOCATE(effpot_mpi%my_rpt)
 end if

end subroutine effpot_mpi_free
!!***

end module m_effpot_mpi
!!***

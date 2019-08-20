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
!! Copyright (C) 2010-2019 ABINIT group (AM)
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
 use m_abicore
 use m_supercell,only: getPBCIndexes_supercell

 implicit none

 public :: effpot_mpi_init
 public :: effpot_mpi_free
!!***

!!****t* m_effpot_mpi/effpot_mpi_type
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
   ! my_cells(4,my_ncell)
   ! 1-3 are the indexes of the cells in the supercell treat by this CPU
   ! dimension 4 is the index of the first atom in the cell

   integer,allocatable :: my_rpt(:)
   ! my_rpt(my_nrpt)
   ! List of rpt treat by this CPU

   integer,allocatable :: my_atmrpt_index(:,:)
   ! my_cells(my_nrpt,my_ncell)
   ! For each cell in the supercell and each rpt, give the index of the first atoms in the rpt cell
   ! Take into acount the PBC.

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
!! index_rpt(3,nrpt) = indexes of the rpt for the ifc
!! sc_size(3) = size of the supercell, 3 3 3 for example
!! natom = natom in the unit cell
!! ndiv  = number of division to consider. For example if ndiv==2,
!!         the mpi will be set over the 2 lvl of parallelisation,
!!         cell and nrpt by considering nrpt / 2 for each CPU
!!         (still experimental, only ndiv=1 available)
!! nrpt  = number of cell to be parallelised for the IFC
!! comm  = MPI communicator
!!
!! OUTPUT
!! effpot_mpi<type(effpot_mpi_type)()> = effpot_mpi datatype
!! PARENTS
!!      m_effective_potential
!!
!! CHILDREN
!!
!! SOURCE

subroutine effpot_mpi_init(index_rpt,sc_size,effpot_mpi,natom,ndiv,nrpt,comm)

  implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm,natom,ndiv,nrpt
!array
 type(effpot_mpi_type), intent(inout) :: effpot_mpi
 integer,intent(in) :: sc_size(3),index_rpt(3,nrpt)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,i3,icell,ii,irpt,irpt_tmp
 integer :: my_rank,ncell_alone,ncell,nproc
 integer :: npcell,nprpt,virt_rank
 integer :: master = 0
 logical :: iam_master = .false.
 character(len=500) :: msg
!array
 integer :: cell_atom2(3)
 integer,allocatable :: rpt_list(:)

! ***********************************************************************

!Set the number of cell in the supercell
 ncell = product(sc_size(:))


 if (any(sc_size <= 0).or.ncell<=0) then
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
 effpot_mpi%my_ncell = int(aint(real(ncell,sp)/(nproc)))

 if(my_rank >= (nproc-ncell_alone)) then
   effpot_mpi%my_ncell = effpot_mpi%my_ncell  + 1
 end if

 if(ndiv>1) effpot_mpi%my_ncell = npcell

!Allocation of array
 ABI_ALLOCATE(effpot_mpi%my_cells,(effpot_mpi%my_ncell))
 ABI_ALLOCATE(effpot_mpi%my_index_cells,(4,effpot_mpi%my_ncell))
 effpot_mpi%my_cells = 0
 effpot_mpi%my_index_cells = 0

 virt_rank = int(aint(real(my_rank,sp)/(ndiv)))

 do icell=1,effpot_mpi%my_ncell
   if(virt_rank >= (nproc-ncell_alone))then
     effpot_mpi%my_cells(icell)=(int(aint(real(ncell,sp)/nproc)))*(virt_rank)+&
&                              (virt_rank - (nproc-ncell_alone)) + icell
   else
     effpot_mpi%my_cells(icell)=(effpot_mpi%my_ncell)*(virt_rank)  + icell
  end if
 end do

 icell = 0
 ii = 0
 do i1 = 1,sc_size(1)
   do i2 = 1,sc_size(2)
     do i3 = 1,sc_size(3)
       ii = ii +1
       if(any(effpot_mpi%my_cells==ii))then
         icell=icell+1
         effpot_mpi%my_index_cells(1,icell) = i1
         effpot_mpi%my_index_cells(2,icell) = i2
         effpot_mpi%my_index_cells(3,icell) = i3
         effpot_mpi%my_index_cells(4,icell) = (effpot_mpi%my_cells(icell)-1)*natom
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
 ABI_ALLOCATE(effpot_mpi%my_atmrpt_index,(effpot_mpi%my_nrpt,effpot_mpi%my_ncell))
 effpot_mpi%my_rpt = 0
 effpot_mpi%my_atmrpt_index = 0


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

 do icell = 1,effpot_mpi%my_ncell
   i1=effpot_mpi%my_index_cells(1,icell)
   i2=effpot_mpi%my_index_cells(2,icell)
   i3=effpot_mpi%my_index_cells(3,icell)
   do irpt_tmp = 1,effpot_mpi%my_nrpt
     irpt = effpot_mpi%my_rpt(irpt_tmp)
!    do irpt = 1,eff_pot%harmonics_terms%ifcs%nrpt
!    get the cell of atom2  (0 0 0, 0 0 1...)
     cell_atom2(1) = i1 + index_rpt(1,irpt)
     cell_atom2(2) = i2 + index_rpt(2,irpt)
     cell_atom2(3) = i3 + index_rpt(3,irpt)
     call getPBCIndexes_supercell(cell_atom2(1:3),sc_size(1:3))
!    index of the second atom in the displacement array
     effpot_mpi%my_atmrpt_index(irpt_tmp,icell) = &
&       ((cell_atom2(1)-1)*sc_size(2)*sc_size(3))*natom+&
&       ((cell_atom2(2)-1)*sc_size(3))*natom+&
&       ((cell_atom2(3)-1))*natom
   end do
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
   effpot_mpi%my_cells(:) = 0
   ABI_DEALLOCATE(effpot_mpi%my_cells)
 end if

 if (allocated(effpot_mpi%my_index_cells)) then
   effpot_mpi%my_index_cells(:,:) = 0
   ABI_DEALLOCATE(effpot_mpi%my_index_cells)
 end if

 if (allocated(effpot_mpi%my_atmrpt_index)) then
   effpot_mpi%my_atmrpt_index(:,:) = 0
   ABI_DEALLOCATE(effpot_mpi%my_atmrpt_index)
 end if

 if (allocated(effpot_mpi%my_rpt)) then
   effpot_mpi%my_rpt(:) = 0
   ABI_DEALLOCATE(effpot_mpi%my_rpt)
 end if

end subroutine effpot_mpi_free
!!***

end module m_effpot_mpi
!!***

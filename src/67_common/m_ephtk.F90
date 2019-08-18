!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_ephtk
!! NAME
!!  m_ephtk
!!
!! FUNCTION
!!  Helper functions common to e-ph calculations.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2018 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_ephtk

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
 use m_dtset

 use m_fstrings,     only : itoa, sjoin, ltoa

 implicit none

 private

 public :: ephtk_set_phmodes_ship     ! Setup a mask to skip accumulating the contribution of certain phonon modes.
 public :: ephtk_set_pertables        ! Set tables for parallelism over perturbations from my_npert and comm
!!***

contains  !=====================================================
!!***

!!****f* m_ephtk/ephtk_set_phmodes_ship
!! NAME
!!  ephtk_set_phmodes_ship
!!
!! FUNCTION
 ! Setup a mask to skip accumulating the contribution of certain phonon modes.
!!
!! INPUT
!!  dtset<dataset_type>=All input variables for this dataset.
!!
!! OUTPUT
!!   phmodes_skip(natom3)) For each mode: 1 to skip this contribution else 0
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine ephtk_set_phmodes_ship(dtset, phmodes_skip)

!Arguments ------------------------------------
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,allocatable,intent(out) :: phmodes_skip(:)

!Local variables ------------------------------
!scalars
 integer :: natom3

! *************************************************************************

 ! Setup a mask to skip accumulating the contribution of certain phonon modes.
 ! By default do not skip, if set skip all but specified
 natom3 = dtset%natom * 3
 ABI_MALLOC(phmodes_skip, (natom3))
 phmodes_skip = 0

 if (all(dtset%eph_phrange /= 0)) then
   if (minval(dtset%eph_phrange) < 1 .or. &
       maxval(dtset%eph_phrange) > natom3 .or. &
       dtset%eph_phrange(2) < dtset%eph_phrange(1)) then
     MSG_ERROR('Invalid range for eph_phrange. Should be between [1, 3*natom] and eph_modes(2) > eph_modes(1)')
   end if
   call wrtout(std_out, sjoin(" Including phonon modes between [", &
               itoa(dtset%eph_phrange(1)), ',', itoa(dtset%eph_phrange(2)), "]"))
   phmodes_skip = 0
   phmodes_skip(dtset%eph_phrange(1):dtset%eph_phrange(2)) = 1
 end if

end subroutine ephtk_set_phmodes_ship
!!***

!!****f* m_ephtk/ephtk_set_pertables
!! NAME
!!  ephtk_set_pertables
!!
!! FUNCTION
!!  Set tables for parallelism over perturbations from my_npert and comm
!!
!! INPUT
!!  natom: Number of atoms
!!  my_npert: Number of atomic perturbations or phonon modes treated by this MPI rank.
!!  comm: MPI communicator for parallelism over atomic perturbations.
!!
!! OUTPUT
!!  integer,allocatable :: my_pinfo(:,:)
!!     my_pinfo(3, my_npert)
!!     my_pinfo(1, ip) gives the `idir` index of the ip-th perturbation.
!!     my_pinfo(2, ip) gives the `ipert` index of the ip-th perturbation.
!!     my_pinfo(3, ip) gives `pertcase`=idir + (ipert-1)*3
!!  integer,allocatable :: pert_table(:,:)
!!     pert_table(2, natom3)
!!     pert_table(1, npert): rank of the processor treating this atomic perturbation.
!!     pert_table(2, npert): imyp index in my_pinfo table, -1 if this rank is not treating ipert.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine ephtk_set_pertables(natom, my_npert, pert_table, my_pinfo, comm)

!Arguments ------------------------------------
 integer,intent(in) :: natom, my_npert, comm
!arrays
 integer,allocatable :: pert_table(:,:)
 integer,allocatable :: my_pinfo(:,:)

!Local variables ------------------------------
!scalars
 integer :: iatom, idir, pertcase, bstart, bstop, ii, ip, natom3, my_rank, nproc
!arrays
 integer :: all_pinfo(3, natom*3)

! *************************************************************************

 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)

 ! Build table with list of perturbations treated by this CPU.
 natom3 = natom * 3
 ABI_MALLOC(my_pinfo, (3, my_npert))
 ABI_MALLOC(pert_table, (2, natom3))

 do iatom=1,natom
   do idir=1,3
     pertcase = idir + (iatom-1) * 3
     all_pinfo(:, pertcase) = [idir, iatom, pertcase]
     pert_table(1, pertcase) = (pertcase - 1) / (natom3 / nproc)
   end do
 end do
 bstart = (natom3 / nproc) * my_rank + 1
 bstop = bstart + my_npert - 1
 my_pinfo = all_pinfo(:, bstart:bstop)

 pert_table(2, :) = -1
 do ii=1,my_npert
   ip = my_pinfo(3, ii)
   pert_table(2, ip) = ii
 end do
 !write(std_out,*)"my_npert", my_npert, "nproc", nproc; write(std_out,*)"my_pinfo", my_pinfo

end subroutine ephtk_set_pertables
!!***

end module m_ephtk
!!***

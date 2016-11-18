!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_kpts
!! NAME
!!  m_kpts
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2008-2016 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
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

MODULE m_kpts

 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_crystal

 implicit none

 private

 public :: kpts_ibz_from_kptrlatt      ! Determines the IBZ, the weights and the BZ from kptrlatt
!!***

!----------------------------------------------------------------------

contains  !============================================================
!!***

!!****f* m_kpts/kpts_ibz_from_kptrlatt
!! NAME
!!  kpts_ibz_from_kptrlatt
!!
!! FUNCTION
!!   Determines the irreducible wedge, the corresponding weights and the list
!!   of k-points in the Brillouin Zone starting from kptrlatt and the set shifts.
!!
!! INPUTS
!!  cryst<crystal_t> = crystalline structure with info on symmetries and time-reversal.
!!  kptrlatt(3,3)=integer coordinates of the primitive vectors of the
!!   lattice reciprocal to the k point lattice to be generated here
!!   If diagonal, the three values are the Monkhorst-Pack usual values, in case of simple cubic.
!!  nshiftk= number of shift vectors in the repeated cell
!!  shiftk(3,nshiftk) = vectors that will be used to determine the shifts from (0. 0. 0.).
!!  [timrev]: Time-reversal symmetry flag. By default this info is passed via crystal.
!!    One can change the default behaviour by passing this option explicitly.
!!    In this case:
!!        if 1, the time reversal operation has to be taken into account.
!!        if 0, no time-reversal symmetry.
!!
!! OUTPUT
!!  nkibz,nkbz = Number of points in IBZ and BZ, respectively.
!!  These arrays are allocated and returned by the routine:
!!  kibz(3,nkibz)=k-points in the IBZ.
!!  wtk(nkibz)=weights of the k-points in the IBZ (normalized to one).
!!  kbz(3,nkbz)=k-points in the BZ.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine kpts_ibz_from_kptrlatt(cryst, kptrlatt, nshiftk, shiftk, nkibz, kibz, wtk, nkbz, kbz, timrev)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'kpts_ibz_from_kptrlatt'
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nshiftk
 integer,intent(out) :: nkibz,nkbz
 integer,optional,intent(in) :: timrev
 type(crystal_t),intent(in) :: cryst
!arrays
 integer,intent(in) :: kptrlatt(3,3)
 real(dp),intent(in) :: shiftk(3,nshiftk)
 real(dp),allocatable,intent(out) :: wtk(:),kibz(:,:),kbz(:,:)

!Local variables-------------------------------
!scalars
 integer,parameter :: option1=1,brav1=1,iout0=0,chksymbreak0=0
 integer :: ik_ibz,max_nkpt,ierr,my_timerev
!arrays
 integer,allocatable :: ibz2bz(:)
 real(dp),allocatable :: wtk_folded(:),wtk_bz(:),my_kbz(:,:)

! *********************************************************************

 ! Call smpbz to get the full grid of k-points `my_kbz`
 ! brav1=1 is able to treat all bravais lattices (same option used in getkgrid)
 max_nkpt = kptrlatt(1,1)*kptrlatt(2,2)*kptrlatt(3,3) &
   + kptrlatt(1,2)*kptrlatt(2,3)*kptrlatt(3,1) &
   + kptrlatt(1,3)*kptrlatt(2,1)*kptrlatt(3,2) &
   - kptrlatt(1,2)*kptrlatt(2,1)*kptrlatt(3,3) &
   - kptrlatt(1,3)*kptrlatt(2,2)*kptrlatt(3,1) &
   - kptrlatt(1,1)*kptrlatt(2,3)*kptrlatt(3,2)
 max_nkpt = max_nkpt * nshiftk

 ABI_MALLOC(my_kbz, (3, max_nkpt))
 call smpbz(brav1, -1, kptrlatt, max_nkpt, nkbz, nshiftk, option1, shiftk, my_kbz)

 ! Call symkpt to get IBZ and weights.
 ABI_MALLOC(ibz2bz, (nkbz))
 ABI_MALLOC(wtk_folded, (nkbz))
 ABI_MALLOC(wtk_bz, (nkbz))
 wtk_bz = one / nkbz ! weights normalized to unity

  ! TODO BE CAREFUL here, the convention used in cryst%timrev is different.
 my_timerev = 1; if (cryst%timrev == 1) my_timerev = 0
 if (present(timrev)) my_timerev = timrev
 call symkpt(chksymbreak0, cryst%gmet, ibz2bz, iout0, my_kbz, nkbz, nkibz, &
   cryst%nsym, cryst%symrec, my_timerev, wtk_bz, wtk_folded)

 ! Allocate output arrays.
 ABI_MALLOC(wtk, (nkibz))
 ABI_MALLOC(kibz, (3, nkibz))
 ABI_MALLOC(kbz, (3, nkbz))
 kbz(:,:) = my_kbz(:, 1:nkbz)

 do ik_ibz=1,nkibz
   wtk(ik_ibz) = wtk_folded(ibz2bz(ik_ibz))
   kibz(:,ik_ibz) = my_kbz(:,ibz2bz(ik_ibz))
 end do
 !ABI_CHECK(abs(sum(wtk) - one) < tol10, "sum(wtk) != one")

 ABI_FREE(ibz2bz)
 ABI_FREE(wtk_folded)
 ABI_FREE(wtk_bz)
 ABI_FREE(my_kbz)

end subroutine kpts_ibz_from_kptrlatt
!!***

!----------------------------------------------------------------------

end module m_kpts
!!***

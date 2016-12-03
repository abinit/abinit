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

 public :: kpts_timrev_from_kptopt     ! Returns the value of timrev from kptopt
 public :: kpts_ibz_from_kptrlatt      ! Determines the IBZ, the weights and the BZ from kptrlatt
!!***

!----------------------------------------------------------------------

contains  !============================================================
!!***

!!****f* m_kpts/kpts_timrev_from_kptopt
!! NAME
!!  kpts_timrev_from_kptopt
!!
!! FUNCTION
!!  Returns the value of timrev from kptopt
!!  1 if the use of time-reversal is allowed; 0 otherwise
!!
!! INPUTS
!!  kptopt=option for the generation of k points
!!    (defines whether spatical symmetries and/or time-reversal can be used)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer pure function kpts_timrev_from_kptopt(kptopt) result(timrev)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'kpts_timrev_from_kptopt'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: kptopt

! *********************************************************************

 timrev = 1; if (any(kptopt == [3, 4])) timrev = 0

end function kpts_timrev_from_kptopt
!!***

!!****f* m_kpts/kpts_ibz_from_kptrlatt
!! NAME
!!  kpts_ibz_from_kptrlatt
!!
!! FUNCTION
!!  Determines the irreducible wedge, the corresponding weights and the list
!!  of k-points in the Brillouin Zone starting from kptrlatt and the set shifts.
!!
!! INPUTS
!!  cryst<crystal_t> = crystalline structure with info on symmetries and time-reversal.
!!  kptopt=option for the generation of k points
!!    (defines whether spatical symmetries and/or time-reversal can be used)
!!  kptrlatt(3,3)=integer coordinates of the primitive vectors of the
!!   lattice reciprocal to the k point lattice to be generated here
!!   If diagonal, the three values are the Monkhorst-Pack usual values, in case of simple cubic.
!!  nshiftk= number of shift vectors in the repeated cell
!!  shiftk(3,nshiftk) = vectors that will be used to determine the shifts from (0. 0. 0.).
!!
!! OUTPUT
!!  nkibz,nkbz = Number of points in IBZ and BZ, respectively.
!!  The following arrays are allocated and returned by the routine:
!!  kibz(3,nkibz) = k-points in the IBZ.
!!  wtk(nkibz) = weights of the k-points in the IBZ (normalized to one).
!!  kbz(3,nkbz) = k-points in the BZ.
!!  [new_kptrlatt] = New value of kptrlatt returned by getkgrid
!!  [new_shiftk(3,new_nshiftk)] = New set of shifts returned by getkgrid
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine kpts_ibz_from_kptrlatt( &
  cryst, kptrlatt, kptopt, nshiftk, shiftk, nkibz, kibz, wtk, nkbz, kbz, &
  new_kptrlatt, new_shiftk)  ! Optional


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'kpts_ibz_from_kptrlatt'
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nshiftk,kptopt
 integer,intent(out) :: nkibz,nkbz
 type(crystal_t),intent(in) :: cryst
!arrays
 integer,intent(in) :: kptrlatt(3,3)
 integer,optional,intent(out) :: new_kptrlatt(3,3)
 real(dp),intent(in) :: shiftk(3,nshiftk)
 real(dp),allocatable,intent(out) :: wtk(:),kibz(:,:),kbz(:,:)
 real(dp),optional,allocatable,intent(out) :: new_shiftk(:,:)

!Local variables-------------------------------
!scalars
 integer,parameter :: iout0=0,chksymbreak0=0,iscf2=2
 integer :: my_nshiftk,nkpt_computed
 real(dp) :: kptrlen
!arrays
 integer,parameter :: vacuum0(3)=[0,0,0]
 integer :: my_kptrlatt(3,3)
 real(dp) :: my_shiftk(3,210)

! *********************************************************************

 ! First call to getkgrid to obtain the number of points in the BZ.
 ABI_MALLOC(kibz, (3,0))
 ABI_MALLOC(wtk, (0))

 ! Copy kptrlatt and shifts because getkgrid can change them
 ! Be careful as getkgrid expects shiftk(3,210).
 ABI_CHECK(nshiftk > 0 .and. nshiftk <= 210, "nshiftk must be in [1,210]")
 my_nshiftk = nshiftk; my_shiftk = zero; my_shiftk(:,1:nshiftk) = shiftk
 my_kptrlatt = kptrlatt

 call getkgrid(chksymbreak0,iout0,iscf2,kibz,kptopt,my_kptrlatt,kptrlen,&
   cryst%nsym,0,nkibz,my_nshiftk,cryst%nsym,cryst%rprimd,my_shiftk,cryst%symafm,cryst%symrel,vacuum0,wtk)

 ABI_FREE(kibz)
 ABI_FREE(wtk)

 ! Recall getkgrid to get kibz and wtk.
 ABI_MALLOC(kibz, (3, nkibz))
 ABI_MALLOC(wtk, (nkibz))

 call getkgrid(chksymbreak0,iout0,iscf2,kibz,kptopt,my_kptrlatt,kptrlen,&
   cryst%nsym,nkibz,nkpt_computed,my_nshiftk,cryst%nsym,cryst%rprimd,my_shiftk,&
   cryst%symafm,cryst%symrel,vacuum0,wtk,fullbz=kbz)

 nkbz = size(kbz, dim=2)

 ! Optionally, return new shifts and new_kptrlatt
 if (present(new_shiftk)) then
   ABI_MALLOC(new_shiftk, (3, my_nshiftk))
   new_shiftk = my_shiftk(:, 1:my_nshiftk)
 end if
 if (present(new_kptrlatt)) new_kptrlatt = my_kptrlatt

 DBG_CHECK(abs(sum(wtk) - one) < tol10, "sum(wtk) != one")

end subroutine kpts_ibz_from_kptrlatt
!!***

!----------------------------------------------------------------------

end module m_kpts
!!***

!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_lgroup
!! NAME
!! m_lgroup
!!
!! FUNCTION
!! The little group of a q-point is defined as the subset of the space group that preserves q,
!! modulo a G0 vector (also called umklapp vector). Namely:
!!
!!    Sq = q + G0
!!
!! where S is an operation in reciprocal space (symrec)
!! If time reversal symmetry can be used, it is possible to enlarge the little group by
!! including the operations such as:
!!
!!   -Sq = q + G0
!!
!! The operations of little group define an irriducible wedge, IBZ(q), that is usually larger
!! than the irredubile zone defined by the point group of the crystal. The two zones coincide when q=0
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2018 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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

module m_lgroup

 use defs_basis
 use m_errors
 use m_abicore
 use m_crystal
 use m_symkpt

 use m_fstrings,   only : ftoa, ktoa, sjoin
 use m_symtk,      only : chkgrp, littlegroup_q

 implicit none

 private
!!***

!!****t* m_sigmaph/lgroup_t
!! NAME
!! lgroup_t
!!
!! FUNCTION
!!  Stores tables associated to the little-group.
!!
!! SOURCE

 type, public :: lgroup_t

   integer :: nibz
   ! Number of points in the IBZ(q)

   integer :: timrev
   ! timrev=1 if time-reversal symmetry can be used, 0 otherwise.

   real(dp) :: point(3)
   ! The external q-point.

   integer,allocatable :: symtab(:,:,:)
   ! symtab(4, 2, cryst%nsym)
   ! nsym is the **total** number of symmetries of the system as given by cryst%nsym
   ! three first numbers define the G vector;
   ! fourth number is zero if the q-vector is not preserved, is 1 otherwise.
   ! second index is one without time-reversal symmetry, two with time-reversal symmetry

   real(dp),allocatable :: ibz(:,:)
   ! ibz(3, nibz)
   ! K-points in the IBZ(q)

   real(dp),allocatable :: weights(:)
   ! weights(nibz)
   ! Weights in the IBZ(q), normalized to 1

 end type lgroup_t
!!***

 public :: lgroup_new      ! Creation method.
 public :: lgroup_print    ! Print the object
 public :: lgroup_free     ! Free memory.

contains  !=====================================================
!!***

!!****f* m_sigmaph/lgroup_new
!! NAME
!!  lgroup_new
!!
!! FUNCTION
!!  Build the little group of the k-point.
!!
!! INPUTS
!!  cryst(crystal_t)=Crystalline structure
!!  kpoint(3)=External k-point defining the little-group
!!  timrev=1 if time-reversal symmetry can be used, 0 otherwise.
!!  nkbz=Number of k-points in the BZ.
!!  kbz(3,nkbz)=K-points in the BZ.
!!  nkibz=Number of k-points in the IBZ
!!  kibz(3,nkibz)=Irreducible zone.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

type (lgroup_t) function lgroup_new(cryst, kpoint, timrev, nkbz, kbz, nkibz, kibz) result(new)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'lgroup_new'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: timrev,nkibz,nkbz
 type(crystal_t),intent(in) :: cryst
!arrays
 real(dp),intent(in) :: kpoint(3),kbz(3,nkbz),kibz(3,nkibz)

!Local variables ------------------------------
!scalars
 integer,parameter :: iout0=0,my_timrev0=0,chksymbreak0=0,debug=0
 integer :: otimrev_k,ierr,itim,isym,nsym_lg,ik
!arrays
 integer :: symrec_lg(3,3,2*cryst%nsym),symafm_lg(3,3,2*cryst%nsym)
 integer,allocatable :: ibz2bz(:)
 real(dp),allocatable :: wtk(:),wtk_folded(:)

! *************************************************************************

 ! TODO: Option to exclude umklapp/time-reversal symmetry and kptopt
 new%point = kpoint
 new%timrev = timrev

 ! Determines the symmetry operations by which the k-point is preserved,
 ABI_MALLOC(new%symtab, (4, 2, cryst%nsym))
 call littlegroup_q(cryst%nsym, kpoint, new%symtab, cryst%symrec, cryst%symafm, otimrev_k, prtvol=0)

 ABI_CHECK(new%timrev==1, "timrev == 0 not coded")
 nsym_lg = 0
 do itim=1,2
   do isym=1,cryst%nsym
     if (cryst%symafm(isym) == -1) cycle
     if (new%symtab(4, itim, isym) /= 1) cycle ! not \pm Sq = q+g0
     nsym_lg = nsym_lg + 1
     symrec_lg(:,:,nsym_lg) = cryst%symrec(:,:,isym) * (-2* itim + 3)
   end do
 end do

 ! Check group closure.
 symafm_lg = 1
 call chkgrp(nsym_lg, symafm_lg, symrec_lg, ierr)
 ABI_CHECK(ierr == 0, "Error in group closure")

 ! Find the irreducible zone with the little group operations.
 ! Do not use time-reversal since it has been manually introduced previously
 ABI_MALLOC(ibz2bz, (nkbz))
 ABI_MALLOC(wtk_folded, (nkbz))
 ABI_MALLOC(wtk, (nkbz))
 wtk = one / nkbz ! Weights sum up to one

 ! TODO: In principle here we would like to have a set that contains the initial IBZ.
 call symkpt(chksymbreak0,cryst%gmet,ibz2bz,iout0,kbz,nkbz,new%nibz,&
   nsym_lg,symrec_lg,my_timrev0,wtk,wtk_folded)

 ABI_MALLOC(new%ibz, (3, new%nibz))
 ABI_MALLOC(new%weights, (new%nibz))

 do ik=1,new%nibz
   new%weights(ik) = wtk_folded(ibz2bz(ik))
   new%ibz(:,ik) = kbz(:, ibz2bz(ik))
 end do
 ABI_CHECK(sum(new%weights) - one < tol12, sjoin("Weights don't sum up to one but to:", ftoa(sum(new%weights))))

 ABI_FREE(ibz2bz)
 ABI_FREE(wtk_folded)
 ABI_FREE(wtk)

 ! Debug section.
 if (debug /= 0) then
   do ik=1,new%nibz
     if (ik <= nkibz) then
       write(std_out,"(a)")sjoin(ktoa(new%ibz(:,ik)), ktoa(new%ibz(:,ik) - kibz(:,ik)), ftoa(new%weights(ik)))
     else
       write(std_out,"(a)")sjoin(ktoa(new%ibz(:,ik)), "[---]", ftoa(new%weights(ik)))
     end if
   end do
 end if

end function lgroup_new
!!***

!----------------------------------------------------------------------

!!****f* m_lgroup/lgroup_print
!! NAME
!! lgroup_print
!!
!! FUNCTION
!!  Print the object
!!
!! INPUTS
!!  [title]=String to be printed as header for additional info.
!!  [unit]=Unit number for output
!!  [prtvol]=Verbosity level
!!
!! OUTPUT
!!  Only printing
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine lgroup_print(self, title, unit, prtvol)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'lgroup_print'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,optional,intent(in) :: unit, prtvol
 character(len=*),optional,intent(in) :: title
 type(lgroup_t),intent(in) :: self

!Local variables-------------------------------
!scalars
 integer :: my_prtvol,my_unt,ik
 character(len=500) :: msg
! *************************************************************************

 my_unt = std_out; if (present(unit)) my_unt = unit
 my_prtvol = 0; if (present(prtvol)) my_prtvol = prtvol

 msg = ' ==== Info on the <lgroup_t> object ==== '
 if (present(title)) msg = ' ==== '//trim(adjustl(title))//' ==== '
 call wrtout(my_unt, msg)

 write(msg, '(3a, 2(a, i0), a)') &
  ' Little group point: ................... ', ktoa(self%point), ch10, &
  ' Number of points in IBZ(p) ............ ', self%nibz, ch10, &
  ' Time-reversal flag (0: No, 1: Yes) .... ', self%timrev
 call wrtout(my_unt, msg)

 if (my_prtvol /= 0) then
   do ik=1,self%nibz
      call wrtout(my_unt, sjoin(ktoa(self%ibz(:,ik)), ftoa(self%weights(ik))))
   end do
 end if

end subroutine lgroup_print
!!***

!!****f* m_sigmaph/lgroup_free
!! NAME
!!  lgroup_free
!!
!! FUNCTION
!!  Free memory
!!
!! PARENTS
!!      m_sigmaph
!!
!! CHILDREN
!!
!! SOURCE

subroutine lgroup_free(self)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'lgroup_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(lgroup_t),intent(inout) :: self

! *************************************************************************

 ! integer
 if (allocated(self%symtab)) then
   ABI_FREE(self%symtab)
 end if
 ! real
 if (allocated(self%ibz)) then
   ABI_FREE(self%ibz)
 end if
 if (allocated(self%weights)) then
   ABI_FREE(self%weights)
 end if

end subroutine lgroup_free
!!***

end module m_lgroup
!!***

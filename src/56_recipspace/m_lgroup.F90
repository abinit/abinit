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
!!  Copyright (C) 2008-2019 ABINIT group (MG)
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
 use m_copy
 use m_symkpt
 use m_sort
 use m_xmpi

 use m_fstrings,      only : ftoa, ktoa, sjoin
 use m_numeric_tools, only : wrap2_pmhalf
 use m_geometry,      only : normv
 use m_kpts,          only : listkk
 use m_symtk,         only : chkgrp, littlegroup_q

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

   integer :: nbz
   ! Number of points in the full BZ

   integer :: nsym_lg
   ! Number of operations in the little group. Including "time-reversed" symmetries
   ! i.e. symmetries S for which -S q = q + G. See also littlegroup_q

   integer :: input_timrev
   ! 1 if time-reversal symmetry can be used, 0 otherwise.
   ! NB: This flag refers to the generation of the initial set of k-points.
   ! One should pay attention when calling other routines in which timrev is required
   ! Because from Sq = q does not necessarily follow that -Sq = q if q is not on zone-border.
   ! The operations in G-space stored here already include time-reversal if input_timrev == 1
   ! so one should call k-point routines with timrev = 0.

   real(dp) :: point(3)
   ! The external q-point.

   integer,allocatable :: symtab(:,:,:)
   ! symtab(4, 2, cryst%nsym)
   ! nsym is the **total** number of spatial symmetries of the system as given by cryst%nsym
   ! three first numbers define the G vector;
   ! fourth number is zero if the q-vector is not preserved, is 1 otherwise.
   ! second index is one without time-reversal symmetry, two with time-reversal symmetry

   integer, allocatable :: symrec_lg(:, :, :)
   ! symrec_lg(3, 3, nsym_lg)
   ! Symmetry operations in G-space (including time-reversed operations if any)

   integer, allocatable :: symafm_lg(:)
   ! symafm_lg(nsym_lg)
   ! Anti-ferromagnetic character

   integer,allocatable :: bz2ibz_smap(:,:)
   ! bz2ibz_smap(nbz, 6) Mapping BZ --> IBZ.
   ! Note that here we used the symmetries of the little group.

   integer, allocatable :: lgsym2glob(:, :)
   ! lgsym2glob(2, nsym_lg)
   ! Mapping isym_lg --> [isym, itime]
   ! where isym is the index of the operaion in crystal%symrec
   ! and itim is 2 if time-reversal T must be included else 1.

  real(dp) :: gmet(3,3)
   ! Reciprocal space metric in bohr^{-2}

   real(dp),allocatable :: ibz(:,:)
   ! ibz(3, nibz)
   ! K-points in the IBZ(q)

   real(dp),allocatable :: weights(:)
   ! weights(nibz)
   ! Weights in the IBZ(q), normalized to 1

 contains

   procedure :: findq_ibzk => lgroup_findq_ibzk
   ! Find the index of the point in the IBZ(k).

   procedure :: find_ibzimage => lgroup_find_ibzimage
   ! Find the symmetrical image in the IBZ(k) of a qpoint in the BZ.

   procedure :: print => lgroup_print
   ! Print the object

   procedure :: free => lgroup_free
   ! Free memory.

 end type lgroup_t
!!***

 public :: lgroup_new                 ! Creation method.

contains  !=====================================================
!!***

!!****f* m_sigmaph/lgroup_new
!! NAME
!!  lgroup_new
!!
!! FUNCTION
!!  Build the little group of the k-point. Return IBZ(k) points packed in shells.
!!  to facilitate optimization of loops.
!!
!! INPUTS
!!  cryst(crystal_t)=Crystalline structure
!!  kpoint(3)=External k-point defining the little-group
!!  timrev=1 if time-reversal symmetry can be used, 0 otherwise.
!!  nkbz=Number of k-points in the BZ.
!!  kbz(3,nkbz)=K-points in the BZ.
!!  nkibz=Number of k-points in the IBZ
!!  kibz(3,nkibz)=Irreducible zone.
!!  comm= MPI communicator.
!!  sord=Defines how to order the points in %ibz. ">" for increasing norm. "<" decreasing. Default: ">"
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

type(lgroup_t) function lgroup_new(cryst, kpoint, timrev, nkbz, kbz, nkibz, kibz, comm, sord) result(new)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: timrev,nkibz,nkbz,comm
 type(crystal_t),intent(in) :: cryst
 character(len=1),optional,intent(in) :: sord
!arrays
 real(dp),intent(in) :: kpoint(3),kbz(3,nkbz),kibz(3,nkibz)

!Local variables ------------------------------
!scalars
 integer,parameter :: iout0=0,my_timrev0=0,chksymbreak0=0,debug=0
 integer :: otimrev_k,ierr,itim,isym,ik_ibz,ik_bz,ksign,isym_lgk
!arrays
 integer :: symrec_lg(3,3,2*cryst%nsym), symafm_lg(2*cryst%nsym), lgsym2glob(2, 2*cryst%nsym)
 real(dp) :: kred(3),shift(3)
 integer,allocatable :: ibz2bz(:), iperm(:), inv_iperm(:)
 real(dp),allocatable :: wtk_folded(:), kord(:,:)

! *************************************************************************

 ! TODO: Option to exclude umklapp/time-reversal symmetry and kptopt
 new%point = kpoint
 new%input_timrev = timrev
 new%gmet = cryst%gmet
 new%nbz = nkbz

 ! Determines the symmetry operations by which the k-point is preserved,
 ABI_MALLOC(new%symtab, (4, 2, cryst%nsym))
 call littlegroup_q(cryst%nsym, kpoint, new%symtab, cryst%symrec, cryst%symafm, otimrev_k, prtvol=0)

 new%nsym_lg = 0
 do itim=1,new%input_timrev + 1
   do isym=1,cryst%nsym
     if (cryst%symafm(isym) == -1) cycle
     if (new%symtab(4, itim, isym) /= 1) cycle ! not \pm Sq = q+g0
     new%nsym_lg = new%nsym_lg + 1
     symrec_lg(:,:,new%nsym_lg) = cryst%symrec(:,:,isym) * (-2* itim + 3)
     symafm_lg(new%nsym_lg) = cryst%symafm(isym)
     lgsym2glob(:, new%nsym_lg) = [isym, itim]
   end do
 end do

 call alloc_copy(symrec_lg(:, :, 1:new%nsym_lg), new%symrec_lg)
 call alloc_copy(symafm_lg(1:new%nsym_lg), new%symafm_lg)
 call alloc_copy(lgsym2glob(:, 1:new%nsym_lg), new%lgsym2glob)

 ! Check group closure.
 if (debug /= 0) then
   call chkgrp(new%nsym_lg, symafm_lg, symrec_lg, ierr)
   ABI_CHECK(ierr == 0, "Error in group closure")
 end if

 ! Find the irreducible zone with the little group operations.
 ! Do not use time-reversal since it has been manually introduced previously
 ABI_MALLOC(ibz2bz, (nkbz))

 ABI_MALLOC(new%bz2ibz_smap, (6, nkbz))
 ! IBZ2BZ ?

 ! TODO: In principle here we would like to have a set that contains the initial IBZ.
 call symkpt_new(chksymbreak0, cryst%gmet, ibz2bz, iout0, kbz, nkbz, new%nibz,&
   new%nsym_lg, new%symrec_lg, my_timrev0, new%bz2ibz_smap, comm)

 ABI_MALLOC(new%ibz, (3, new%nibz))
 ABI_CALLOC(new%weights, (new%nibz))

 do ik_bz=1,nkbz
   ik_ibz   = new%bz2ibz_smap(1,ik_bz)
   isym_lgk = new%bz2ibz_smap(2,ik_bz)
   new%bz2ibz_smap(2,ik_bz) = lgsym2glob(1,isym_lgk)
   new%bz2ibz_smap(3,ik_bz) = lgsym2glob(2,isym_lgk)
   new%weights(ik_ibz) = new%weights(ik_ibz) + 1
 end do
 new%weights(:) = new%weights(:) / nkbz

 do ik_ibz=1,new%nibz
   ik_bz = ibz2bz(ik_ibz)
   new%ibz(:,ik_ibz) = kbz(:, ik_bz)
 end do

 ! TODO: Activate this part so that we can cache the q-point in the IBZ.
 ! Results are ok but this change is postponed because it leads to an increase
 ! in the walltime spent in listkk likely because of the different order.

 ! Need to repack the IBZ points and rearrange the other arrays dimensioned with nibz.
 ! In principle, the best approach would be to pack in stars using crystal%symrec.
 ! For the time being we pack in shells (much easier). Use wtk_folded as workspace to store the norm.
 ksign = 0
 if (present(sord)) then
   if (sord == "<") ksign = -1
   if (sord == ">") ksign = +1
 end if

 if (ksign /= 0) then
   ABI_MALLOC(wtk_folded, (new%nibz))
   do ik_ibz=1,new%nibz
     call wrap2_pmhalf(new%ibz(:, ik_ibz), kred, shift)
     wtk_folded(ik_ibz) = ksign * normv(kred, cryst%gmet, "G")
   end do

   ABI_MALLOC(iperm, (new%nibz))
   iperm = [(ik_ibz, ik_ibz=1, new%nibz)]
   call sort_dp(new%nibz, wtk_folded, iperm, tol12)
   !iperm = [(ik_ibz, ik_ibz=1, new%nibz)]

   ! Trasfer data.
   ABI_MALLOC(kord, (3, new%nibz))
   do ik_ibz=1,new%nibz
     kord(:, ik_ibz) = new%ibz(:, iperm(ik_ibz))
     wtk_folded(ik_ibz) = new%weights(iperm(ik_ibz))
   end do
   new%ibz = kord(:, 1:new%nibz)
   new%weights = wtk_folded(1:new%nibz)

   ! Rearrange bz2ibz_smap as well --> need the inverse of iperm.
   ABI_MALLOC(inv_iperm, (new%nibz))
   do ik_ibz=1,new%nibz
     inv_iperm(iperm(ik_ibz)) = ik_ibz
   end do

   do ik_bz=1,new%nbz
     ik_ibz = new%bz2ibz_smap(1, ik_bz)
     new%bz2ibz_smap(1, ik_bz) = inv_iperm(ik_ibz)
   end do

   ABI_FREE(inv_iperm)
   ABI_FREE(kord)
   ABI_FREE(iperm)
   ABI_FREE(wtk_folded)
 end if

 ABI_FREE(ibz2bz)

 ! Debug section.
 ABI_CHECK(sum(new%weights) - one < tol6, sjoin("Weights don't sum up to one but to:", ftoa(sum(new%weights))))

 if (debug /= 0) then
   do ik_ibz=1,new%nibz
     if (ik_ibz <= nkibz) then
       write(std_out,"(a)")sjoin(ktoa(new%ibz(:,ik_ibz)), ktoa(new%ibz(:,ik_ibz) - kibz(:,ik_ibz)), ftoa(new%weights(ik_ibz)))
     else
       write(std_out,"(a)")sjoin(ktoa(new%ibz(:,ik_ibz)), "[---]", ftoa(new%weights(ik_ibz)))
     end if
   end do
 end if

end function lgroup_new
!!***

!!****f* m_lgroup/lgroup_findq_ibzk
!! NAME
!!  lgroup_findq_ibzk
!!
!! FUNCTION
!!  Find the index of the q-point in the IBZ(k). Umklapp vectors are not allowed.
!!  Returns -1 if not found.
!!
!! INPUTS
!!  qpt(3)=q-point in reduced coordinates.
!!  [qtol]=Optional tolerance for q-point comparison.
!!         For each reduced direction the absolute difference between the coordinates must be less that qtol
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer pure function lgroup_findq_ibzk(self, qpt, qtol) result(iqpt)

!Arguments ------------------------------------
!scalars
 real(dp),optional,intent(in) :: qtol
 class(lgroup_t),intent(in) :: self
!arrays
 real(dp),intent(in) :: qpt(3)

!Local variables-------------------------------
 integer :: iq
 real(dp) :: my_qtol

! *************************************************************************

 my_qtol = tol6; if (present(qtol)) my_qtol = qtol

 iqpt = -1
 do iq=1,self%nibz
   if (all(abs(self%ibz(:, iq) - qpt) < my_qtol)) then
      iqpt = iq; exit
   end if
 end do

end function lgroup_findq_ibzk
!!***

!----------------------------------------------------------------------

!!****f* m_lgroup/lgroup_find_ibzimage
!! NAME
!! lgroup_find_ibzimage
!!
!! FUNCTION
!!  Find the symmetrical image in the IBZ(k) of a qpoint in the BZ
!!  Returns -1 if not found.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function lgroup_find_ibzimage(self, qpt) result(iq_ibz)

!Arguments ------------------------------------
 class(lgroup_t),intent(in) :: self
 real(dp),intent(in) :: qpt(3)

!Local variables-------------------------------
!scalars
 integer, parameter :: timrev0 = 0
 real(dp) :: dksqmax
!arrays
 integer :: indkk(6)
! *************************************************************************

 ! Note use_symrec and timrev0
 call listkk(dksqmax, self%gmet, indkk, self%ibz, qpt, self%nibz, 1, self%nsym_lg, &
    1, self%symafm_lg, self%symrec_lg, timrev0, xmpi_comm_self, exit_loop=.True., use_symrec=.True.)

 iq_ibz = indkk(1)
 if (dksqmax > tol12) iq_ibz = -1

end function lgroup_find_ibzimage
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

!Arguments ------------------------------------
 integer,optional,intent(in) :: unit, prtvol
 character(len=*),optional,intent(in) :: title
 class(lgroup_t),intent(in) :: self

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

 write(msg, '(3a, 2(a, i0, a))') &
  ' Little group point: ................... ', trim(ktoa(self%point)), ch10, &
  ' Number of points in IBZ(p) ............ ', self%nibz, ch10, &
  ' Time-reversal flag (0: No, 1: Yes) .... ', self%input_timrev, ch10
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

!Arguments ------------------------------------
 class(lgroup_t),intent(inout) :: self

! *************************************************************************

 ! integer
 ABI_SFREE(self%symrec_lg)
 ABI_SFREE(self%symafm_lg)
 ABI_SFREE(self%bz2ibz_smap)
 ABI_SFREE(self%lgsym2glob)
 ABI_SFREE(self%symtab)
 ! real
 ABI_SFREE(self%ibz)
 ABI_SFREE(self%weights)

end subroutine lgroup_free
!!***

end module m_lgroup
!!***

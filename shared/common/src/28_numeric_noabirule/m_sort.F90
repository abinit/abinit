!!****m* ABINIT/m_sort
!! NAME
!! m_sort
!!
!! FUNCTION
!! Sorting algorithms.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2026 ABINIT group (XG, MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_sort

 use defs_basis
 use m_errors
 use m_profiling_abi

 implicit none

 private

 ! Low-level routines
 public :: sort_dp       ! Sort double precision array
 public :: sort_int      ! Sort integer array

 ! Helper functions to perform common operations.
 public :: sort_rpts     ! Sort list of real points by |r|
 public :: sort_rvals    ! Out-of-place sort of real values
 public :: sort_gvecs    ! Sort list of g-vectors by norm

CONTAINS  !====================================================================================================
!!***

!!****f* m_sort/sort_dp
!! NAME
!!  sort_dp
!!
!! FUNCTION
!!  Sort double precision array list(n) into ascending numerical order using Heapsort
!!  algorithm, while making corresponding rearrangement of the integer
!!  array iperm. Consider that two double precision numbers within tolerance tol are equal.
!!
!! INPUTS
!!  n: dimension of the list
!!  tol: numbers within tolerance are equal
!!  list(n)  intent(inout) list of double precision numbers to be sorted
!!  iperm(n) intent(inout) iperm(i)=i (very important)
!!  [order]: order of sorting (ascending or descending)
!!
!! OUTPUT
!!  list(n)  sorted list
!!  iperm(n) index of permutation giving the right ascending order:
!!      the i-th element of the ouput ordered list had index iperm(i) in the input list.
!!
!! SOURCE

subroutine sort_dp(n, list, iperm, tol, order)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: n
 integer, intent(inout) :: iperm(n)
 real(dp), intent(inout) :: list(n)
 real(dp), intent(in) :: tol
 integer, optional, intent(in) :: order

!Local variables-------------------------------
!scalars
 integer :: l,ir,iap,i,j,my_order
 real(dp) :: ap
 character(len=500) :: msg

 if (n==1) then

! Accomodate case of array of length 1: already sorted!
  return

 else if (n<1) then

! Should not call with n<1
  write(msg, "(a,i0,2a)")&
    "sort_dp has been called with array length n= ",n, ch10, &
    "having a value less than 1. This is not allowed."
  ABI_ERROR(msg)

 else ! n>1

! Conduct the usual sort

  l=n/2+1
  ir=n

  my_order = 1; if (present(order)) my_order = order
  if (abs(my_order) /= 1) then
    write(msg, "(a,i0,2a)")&
      "sort_dp has been called with an invalid order= ",my_order, ch10, &
      "This is not allowed."
    ABI_ERROR(msg)
  end if
  if (my_order==-1) then
    ! Descending order
    list = -list
  end if

  do   ! Infinite do-loop

   if (l>1) then

    l=l-1
    ap=list(l)
    iap=iperm(l)

   else ! l<=1

    ap=list(ir)
    iap=iperm(ir)
    list(ir)=list(1)
    iperm(ir)=iperm(1)
    ir=ir-1

    if (ir==1) then
     list(1)=ap
     iperm(1)=iap
     exit   ! This is the end of this algorithm
    end if

   end if ! l>1

   i=l
   j=l+l

   do while (j<=ir)
    if (j<ir) then
     if ( list(j)<list(j+1)-tol .or.  &
&        (list(j)<list(j+1)+tol.and.iperm(j)<iperm(j+1))) j=j+1
    endif
    if (ap<list(j)-tol .or. (ap<list(j)+tol.and.iap<iperm(j))) then
     list(i)=list(j)
     iperm(i)=iperm(j)
     i=j
     j=j+j
    else
     j=ir+1
    end if
   enddo

   list(i)=ap
   iperm(i)=iap

  enddo ! End infinite do-loop

  if (my_order == -1) then
    ! Descending order
    list = -list
  end if

 end if ! n>1

end subroutine sort_dp
!!***

!!****f* m_sort/sort_int
!! NAME
!!  sort_int
!!
!! FUNCTION
!!   Sort integer array list(n) into ascending numerical order using Heapsort
!!   algorithm, while making corresponding rearrangement of the integer array iperm.
!!
!! INPUTS
!!  n: dimension of the list
!!  list(n)  intent(inout) list of double precision numbers to be sorted
!!  iperm(n) intent(inout) iperm(i)=i (very important)
!!  [order]: order of sorting (ascending or descending)
!!
!! OUTPUT
!!  list(n): sorted list
!!  iperm(n): index of permutation given the right ascending order
!!      the i-th element of the ouput ordered list had index iperm(i) in the input list.
!!
!! SOURCE

subroutine sort_int(n, list, iperm, order)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
 integer,intent(inout) :: list(n),iperm(n)
 integer,optional,intent(in) :: order

!Local variables-------------------------------
!scalars
 integer :: l,ir,i,j,ip,ipp,my_order
 character(len=500) :: msg
! *************************************************************************

 if (n==1) then

! Accomodate case of array of length 1: already sorted!
  return

 else if (n<1) then

! Should not call with n<1
  write(msg, "(a,i0,2a)")&
    "sort_int has been called with array length n= ",n, ch10, &
    "having a value less than 1. This is not allowed."
  ABI_ERROR(msg)

 else ! n>1

! Conduct the usual sort

  l=n/2+1
  ir=n

  my_order = 1; if (present(order)) my_order = order
  if (abs(my_order) /= 1) then
    write(msg, "(a,i0,2a)")&
      "sort_int has been called with an invalid order= ",my_order, ch10, &
      "This is not allowed."
    ABI_ERROR(msg)
  end if
  if (my_order==-1) then
    ! Descending order
    list = -list
  end if

  do   ! Infinite do-loop

   if (l>1) then

    l=l-1
    ip=list(l)
    ipp=iperm(l)

   else

    ip=list(ir)
    ipp=iperm(ir)
    list(ir)=list(1)
    iperm(ir)=iperm(1)
    ir=ir-1

    if (ir==1) then
     list(1)=ip
     iperm(1)=ipp
     exit   ! This is the end of this algorithm
    end if

   end if ! l>1

   i=l
   j=l+l

   do while (j<=ir)
    if (j<ir) then
     if (list(j).lt.list(j+1)) j=j+1
    end if
    if (ip.lt.list(j)) then
     list(i)=list(j)
     iperm(i)=iperm(j)
     i=j
     j=j+j
    else
     j=ir+1
    end if
   enddo

   list(i)=ip
   iperm(i)=ipp

  enddo ! End infinite do-loop

  if (my_order == -1) then
    ! Descending order
    list = -list
  end if

 end if ! n>1

end subroutine sort_int
!!***

!----------------------------------------------------------------------

!!****f* m_sort/sort_rpts
!! NAME
!!  sort_rpts
!!
!! FUNCTION
!!  Sort list of real space 3d-points by norm (ascending order)
!!  Input list is not modified.
!!
!! INPUTS
!!  n: dimension of the list
!!  rpts(3, n): points in reduced coordinates
!!  metric: Metric used to compute |r|.
!!  [tol]: numbers within tolerance are equal.
!!  [order]: order of sorting (ascending or descending)
!!
!! OUTPUT
!!  iperm(n) index of permutation giving the right ascending order:
!!      the i-th element of the ordered list had index iperm(i) in rpts.
!!  [rmod(n)]= list of sorted |r| values.
!!
!! SOURCE

subroutine sort_rpts(n, rpts, metric, iperm, tol, rmod, order)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
 integer,allocatable,intent(out) :: iperm(:)
 real(dp),optional,allocatable,intent(out) :: rmod(:)
 real(dp),optional,intent(in) :: tol
 integer,optional,intent(in) :: order
!arrays
 real(dp),intent(in) :: rpts(3,n), metric(3,3)

!Local variables-------------------------------
!scalars
 integer :: ii, my_order
 real(dp) :: my_tol
!arrays
 real(dp),allocatable :: my_rmod(:)

!************************************************************************

 my_tol = tol12; if (present(tol)) my_tol = tol
 my_order = 1; if (present(order)) my_order = order

 ABI_MALLOC(my_rmod, (n))
 do ii=1,n
   my_rmod(ii) = sqrt(dot_product(rpts(:, ii), matmul(metric, rpts(:, ii))))
 end do
 ABI_MALLOC(iperm, (n))
 iperm = [(ii, ii=1,n)]
 call sort_dp(n, my_rmod, iperm, my_tol, order = my_order)

 if (present(rmod)) then
   call move_alloc(my_rmod, rmod)
 else
   ABI_FREE(my_rmod)
 end if

end subroutine sort_rpts
!!***

!----------------------------------------------------------------------

!!****f* m_sort/sort_rvals
!! NAME
!!  sort_rvals
!!
!! FUNCTION
!!  Sort list of real values (ascending order)
!!  Input list is not modified.
!!
!! INPUTS
!!  n: dimension of the list
!!  in_vals(n): input weigts.
!!  [tol]: tolerance for comparison
!!  [order]: order of sorting (ascending or descending)
!!
!! OUTPUT
!!  iperm(n) index of permutation giving the right ascending order:
!!      the i-th element of the ordered list had index iperm(i) in in_vals.
!!  [sorted_in_vals(n)]= list of sorted values.
!!
!! SOURCE

subroutine sort_rvals(n, in_vals, iperm, sorted_vals, tol, order)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
 real(dp),intent(in) :: in_vals(n)
 integer,allocatable,intent(out) :: iperm(:)
 real(dp),allocatable,intent(out) :: sorted_vals(:)
 real(dp),optional,intent(in) :: tol
 integer,optional,intent(in) :: order

!Local variables-------------------------------
!scalars
 integer :: ii, my_order
 real(dp) :: my_tol

!************************************************************************

 my_tol = tol12; if (present(tol)) my_tol = tol
 my_order = 1; if (present(order)) my_order = order

 ABI_MALLOC(sorted_vals, (n))
 sorted_vals = in_vals
 ABI_MALLOC(iperm, (n))
 iperm = [(ii, ii=1,n)]
 call sort_dp(n, sorted_vals, iperm, my_tol, order = my_order)

end subroutine sort_rvals
!!***

!----------------------------------------------------------------------

!!****f* m_sort/sort_gvecs
!! NAME
!!  sort_gvecs
!!
!! FUNCTION
!!  Sort list of g-vectors (ascending order)
!!  Input list is not modified.
!!
!! INPUTS
!!  npw_k: dimension of the list
!!  kpoint(3): K-point
!!  gmet(3,3): metric matrix.
!!  kg_k(3,npw_k): input weigts.
!!  [tol]: tolerance for comparison
!!  [order]: order of sorting (ascending or descending)
!!
!! OUTPUT
!!  [out_gvec(3,npw_k)]: list of sorted g-vectors
!!  [iperm(npw_k) index of permutation giving the right ascending order:
!!      the i-th element of the ordered list had index iperm(i) in in_vals.
!!
!! SOURCE

subroutine sort_gvecs(npw_k, kpoint, gmet, in_gvec, out_gvec, iperm, tol, order)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw_k, in_gvec(3,npw_k)
 real(dp),intent(in) :: kpoint(3), gmet(3,3)
 integer,allocatable,intent(out) :: out_gvec(:,:)
 integer,allocatable,optional,intent(out) :: iperm(:)
 real(dp),optional,intent(in) :: tol
 integer,optional,intent(in) :: order

!Local variables-------------------------------
!scalars
 integer :: ig, ig_sort, my_order
 real(dp) :: my_tol
 integer,allocatable :: iperm__(:)
 real(dp),allocatable :: kin_kg(:)

!************************************************************************

 my_tol = tol14; if (present(tol)) my_tol = tol
 my_order = 1; if (present(order)) my_order = order

 ABI_MALLOC(kin_kg, (npw_k))
 ABI_MALLOC(iperm__, (npw_k))
 iperm__ = [(ig, ig=1,npw_k)]
 do ig=1,npw_k
   kin_kg(ig) = half * normv(kpoint + in_gvec(:, ig), gmet, "G") ** 2
 end do

 call sort_dp(npw_k, kin_kg, iperm__, my_tol, order = my_order)
 ABI_FREE(kin_kg)

 ABI_MALLOC(out_gvec, (3, npw_k))
 do ig=1,npw_k
   ig_sort = iperm__(ig)
   out_gvec(:,ig) = in_gvec(:,ig_sort)
 end do

 if (present(iperm)) then
   ABI_MALLOC(iperm, (npw_k))
   iperm = iperm__
 end if
 ABI_FREE(iperm__)

contains
function normv(xv, met, space) result(res)

!Arguments ------------------------------------
!scalars
 real(dp) :: res
 character(len=1),intent(in) :: space
!arrays
 real(dp),intent(in) :: met(3,3)
 real(dp),intent(in) :: xv(3)

! *************************************************************************

 res =  ( xv(1)*met(1,1)*xv(1) + xv(2)*met(2,2)*xv(2) + xv(3)*met(3,3)*xv(3)  &
&  +two*( xv(1)*met(1,2)*xv(2) + xv(1)*met(1,3)*xv(3) + xv(2)*met(2,3)*xv(3)) )

 select case (space)
 case ('r','R')
   res=SQRT(res)
 case ('g','G')
   res=two_pi*SQRT(res)
 case default
   ABI_BUG('Wrong value for space')
 end select

end function normv
!!***

end subroutine sort_gvecs
!!***

end module m_sort
!!***

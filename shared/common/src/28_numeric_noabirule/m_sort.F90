!!****m* ABINIT/m_sort
!! NAME
!! m_sort
!!
!! FUNCTION
!! Sorting algorithms.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2020 ABINIT group (XG)
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
 !
 public :: sort_rpts     ! Sort list of real points by |r|
 public :: sort_weights  ! Sort list of weights.

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
!!
!! OUTPUT
!!  list(n)  sorted list
!!  iperm(n) index of permutation giving the right ascending order:
!!      the i-th element of the ouput ordered list had index iperm(i) in the input list.
!!
!! PARENTS
!!      m_bader,m_bz_mesh,m_chi0tk,m_cut3d,m_ebands,m_epjdos,m_exc_diago
!!      m_geometry,m_gsphere,m_ifc,m_ingeo,m_io_screening,m_kpts,m_lgroup
!!      m_mkcore,m_mkrho,m_mlwfovlp_qp,m_mpi_setup,m_outscfcv,m_paw_mkrho
!!      m_paw_pwaves_lmn,m_phonons,m_polynomial_coeff,m_screen,m_skw,m_sort
!!      m_symkpt,m_tddft,m_thmeig,m_use_ga,m_vcoul,m_wvl_rho,mkcore_wvl
!!
!! CHILDREN
!!      move_alloc,sort_dp
!!
!! SOURCE

subroutine sort_dp(n, list, iperm, tol)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: n
 integer, intent(inout) :: iperm(n)
 real(dp), intent(inout) :: list(n)
 real(dp), intent(in) :: tol

!Local variables-------------------------------
!scalars
 integer :: l,ir,iap,i,j
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
!!
!! OUTPUT
!!  list(n): sorted list
!!  iperm(n): index of permutation given the right ascending order
!!      the i-th element of the ouput ordered list had index iperm(i) in the input list.
!!
!! PARENTS
!!      m_dvdb,m_elphon,m_fftcore,m_geometry,m_hdr,m_invars2,m_mpi_setup
!!      m_mpinfo,m_nesting,m_rec,m_spacepar,m_wfk
!!
!! CHILDREN
!!      move_alloc,sort_dp
!!
!! SOURCE

subroutine sort_int(n,list,iperm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
 integer,intent(inout) :: list(n),iperm(n)

!Local variables-------------------------------
!scalars
 integer :: l,ir,i,j,ip,ipp
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
!!
!! OUTPUT
!!  iperm(n) index of permutation giving the right ascending order:
!!      the i-th element of the ordered list had index iperm(i) in rpts.
!!  [rmod(n)]= list of sorted |r| values.
!!
!! PARENTS
!!      m_dvdb,m_sigmaph
!!
!! CHILDREN
!!      move_alloc,sort_dp
!!
!! SOURCE

subroutine sort_rpts(n, rpts, metric, iperm, tol, rmod)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
 integer,allocatable,intent(out) :: iperm(:)
 real(dp),optional,allocatable,intent(out) :: rmod(:)
 real(dp),optional,intent(in) :: tol
!arrays
 real(dp),intent(in) :: rpts(3,n), metric(3,3)

!Local variables-------------------------------
!scalars
 integer :: ii
 real(dp) :: my_tol
!arrays
 real(dp),allocatable :: my_rmod(:)

!************************************************************************

 my_tol = tol12; if (present(tol)) my_tol = tol

 ABI_MALLOC(my_rmod, (n))
 do ii=1,n
   my_rmod(ii) = sqrt(dot_product(rpts(:, ii), matmul(metric, rpts(:, ii))))
 end do
 ABI_MALLOC(iperm, (n))
 iperm = [(ii, ii=1,n)]
 call sort_dp(n, my_rmod, iperm, my_tol)

 if (present(rmod)) then
   call move_alloc(my_rmod, rmod)
 else
   ABI_FREE(my_rmod)
 end if

end subroutine sort_rpts
!!***

!----------------------------------------------------------------------

!!****f* m_sort/sort_weights
!! NAME
!!  sort_weights
!!
!! FUNCTION
!!  Sort list of real values (ascending order)
!!  Input list is not modified.
!!
!! INPUTS
!!  n: dimension of the list
!!  weights(n): input weigts.
!!  [tol]: numbers within tolerance are equal.
!!
!! OUTPUT
!!  iperm(n) index of permutation giving the right ascending order:
!!      the i-th element of the ordered list had index iperm(i) in weights.
!!  [sorted_weights(n)]= list of sorted weigts.
!!
!! PARENTS
!!
!! CHILDREN
!!      move_alloc,sort_dp
!!
!! SOURCE

subroutine sort_weights(n, weights, iperm, tol, sorted_weights)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
 integer,allocatable,intent(out) :: iperm(:)
 real(dp),optional,allocatable,intent(out) :: sorted_weights(:)
 real(dp),optional,intent(in) :: tol
!arrays
 real(dp),intent(in) :: weights(n)

!Local variables-------------------------------
!scalars
 integer :: ii
 real(dp) :: my_tol
!arrays
 real(dp),allocatable :: my_weights(:)

!************************************************************************

 my_tol = tol12; if (present(tol)) my_tol = tol

 ABI_MALLOC(my_weights, (n))
 my_weights = weights
 ABI_MALLOC(iperm, (n))
 iperm = [(ii, ii=1,n)]
 call sort_dp(n, my_weights, iperm, my_tol)

 if (present(sorted_weights)) then
   call move_alloc(my_weights, sorted_weights)
 else
   ABI_FREE(my_weights)
 end if

end subroutine sort_weights
!!***

end module m_sort
!!***

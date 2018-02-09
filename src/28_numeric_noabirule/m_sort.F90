!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_sort
!! NAME
!! m_sort
!!
!! FUNCTION
!! Sorting algorithms.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2018 ABINIT group (XG)
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

 use defs_basis, only : std_out
 use m_errors

 implicit none

 private 

 public :: sort_dp       ! Sort double precision array
 public :: sort_int      ! Sort integer array

CONTAINS  !====================================================================================================
!!***

!!****f* m_sort/sort_dp
!! NAME
!!  sort_dp
!!
!! FUNCTION 
!!  Sort double precision array list(n) into ascending numerical order using Heapsort
!!  algorithm, while making corresponding rearrangement of the integer
!!  array iperm. Consider that two double precision numbers
!!  within tolerance tol are equal.
!!
!! INPUTS
!!  n        intent(in)    dimension of the list
!!  tol      intent(in)    numbers within tolerance are equal
!!  list(n)  intent(inout) list of double precision numbers to be sorted
!!  iperm(n) intent(inout) iperm(i)=i (very important)
!!
!! OUTPUT
!!  list(n)  sorted list
!!  iperm(n) index of permutation given the right ascending order
!!
!! PARENTS
!!      atomden,cpdrv,critics,denfgr,finddistrproc,invacuum,listkk,m_bz_mesh
!!      m_chi0,m_cut3d,m_exc_diago,m_gsphere,m_ifc,m_io_screening
!!      m_paw_pwaves_lmn,m_phonons,m_polynomial_coeff,m_screen,m_skw,m_use_ga
!!      m_vcoul,mkcore,mkcore_inner,mkcore_wvl,mlwfovlp_qp,outscfcv
!!      partial_dos_fractions,shellstruct,symkpt,tddft,thmeig,wvl_initro
!!
!! CHILDREN
!!
!! SOURCE

subroutine sort_dp(n,list,iperm,tol)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sort_dp'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: n
 integer, intent(inout) :: iperm(n)
 double precision, intent(inout) :: list(n)
 double precision, intent(in) :: tol

!Local variables-------------------------------
!scalars
 integer :: l,ir,iap,i,j
 double precision :: ap

 if (n==1) then

! Accomodate case of array of length 1: already sorted!
  return

 else if (n<1) then

! Should not call with n<1
  write(std_out,1000) n
  1000  format(/,' sort_dp has been called with array length n=',i12,/, &
&  ' having a value less than 1.  This is not allowed.')
  MSG_ERROR("Aborting now")

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
!!   algorithm, while making corresponding rearrangement of the integer
!!   array iperm. 
!!
!! INPUTS
!!  n        intent(in)    dimension of the list
!!  list(n)  intent(inout) list of double precision numbers to be sorted
!!  iperm(n) intent(inout) iperm(i)=i (very important)
!!
!! OUTPUT
!!  list(n)  sorted list
!!  iperm(n) index of permutation given the right ascending order
!!
!! PARENTS
!!      getng,getngrec,initmpi_img,invars2,irrzg,m_dvdb,m_hdr,m_nesting,m_wfk
!!      mkfskgrid,shellstruct
!!
!! CHILDREN
!!
!! SOURCE

subroutine sort_int(n,list,iperm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sort_int'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n 
 integer,intent(inout) :: list(n),iperm(n)

!Local variables-------------------------------
!scalars
 integer :: l,ir,i,j,ip,ipp
! *************************************************************************

 if (n==1) then

! Accomodate case of array of length 1: already sorted!
  return

 else if (n<1) then

! Should not call with n<1
  write(std_out,1000) n
  1000  format(/,' sort_int has been called with array length n=',i12,/, &
&  ' having a value less than 1.  This is not allowed.')
  MSG_ERROR("Aborting now")

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

end module m_sort
!!***

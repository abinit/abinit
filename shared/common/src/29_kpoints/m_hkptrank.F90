!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_hkptrank
!! NAME
!! m_hkptrank
!!
!! FUNCTION
!! This module deals with rank objects for hashing k-point vector lists
!! HM: This is based on the kptrank implementation by MJV but using an hashtable to
!! store the inverse mapping and save a lot of memory.
!!
!! COPYRIGHT
!! Copyright (C) 2010-2019 ABINIT group (HM,MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! PARENTS
!!
!! SOURCE


#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_hkptrank

 use defs_basis
 use m_errors
 use m_numeric_tools,   only : wrap2_zero_one
 use m_hashtable

 implicit none

 private

!!***

!!****t* m_hkptrank/hkptrank_type
!! NAME
!! hkptrank_type
!!
!! FUNCTION
!!  structure to contain a rank/inverse rank pair of arrays, with dimensions
!!
!! SOURCE

 type,public :: hkptrank_t
   integer :: max_linear_density
   integer :: max_rank
   integer :: nkpt
   type(hashtable_t) :: invrank
 end type hkptrank_t

 public :: hkptrank_init      ! Sets up the kpt ranks for comparing kpts
 public :: hkptrank_get_1kpt  ! Calculates the rank for one kpt
 public :: hkptrank_get_index ! Get the index of a k-point
 public :: hkptrank_free      ! Free memory
!!***

contains
!!***

!!****f* m_hkptrank/hkptrank_init
!!
!! NAME
!! hkptrank_init
!!
!! FUNCTION
!! This routine sets up the kpt ranks for comparing kpts
!!
!! INPUTS
!!  npt = number of kpoints (eventually irreducible)
!!  kpt = coordinates of kpoints
!!
!! OUTPUT
!!  krank = object containing ranking and inverse ranking
!!
!! PARENTS
!!      get_full_kgrid,m_ddk,m_ebands,m_fstab,m_nesting,m_phgamma,m_pptools
!!      m_tetrahedron,mkfskgrid,mkqptequiv,order_fs_kpts,outelph,read_el_veloc
!!
!! CHILDREN
!!
!! SOURCE

subroutine hkptrank_init(krank,kpt,nkpt)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt
!arrays
 type(hkptrank_t), intent(out) :: krank
 double precision,intent(in) :: kpt(3,nkpt)

!Local variables -------------------------
!scalars
 integer :: ikpt, irank
 double precision :: smallestlen

! *********************************************************************

! find smallest linear length
 smallestlen = one
 do ikpt=1, nkpt
   if (abs(kpt(1,ikpt)) > tol10) smallestlen = min(smallestlen, abs(kpt(1,ikpt)))
   if (abs(kpt(2,ikpt)) > tol10) smallestlen = min(smallestlen, abs(kpt(2,ikpt)))
   if (abs(kpt(3,ikpt)) > tol10) smallestlen = min(smallestlen, abs(kpt(3,ikpt)))
 end do

 krank%max_linear_density = int(one/smallestlen)+1
 krank%max_rank = 2*krank%max_linear_density**3
 krank%nkpt = nkpt
 krank%invrank = hashtable_init(nkpt,4,4)

 ! the following fills invrank for the k-points in the list provided (may be only the irred kpts)
 do ikpt=1,nkpt
   irank = hkptrank_get_1kpt(krank,kpt(:,ikpt))
   call hashtable_add(krank%invrank,irank,ikpt)
 end do

 call hashtable_cleanup(krank%invrank)

end subroutine hkptrank_init
!!***

!----------------------------------------------------------------------

!!****f* m_hkptrank/hkptrank_get_1kpt
!!
!! NAME
!! hkptrank_get_1kpt
!!
!! FUNCTION
!! This routine calculates the rank for one kpt
!!
!! INPUTS
!!  kpt = coordinates of kpoints
!!  krank = rank object for the k-grid we are using
!!
!! OUTPUT
!!  rank = rank of the kpoint
!!
!! PARENTS
!!      elphon,get_full_kgrid,integrate_gamma,integrate_gamma_alt,k_neighbors
!!      m_ddk,m_kptrank,m_nesting,m_pptools,m_tetrahedron,mkfskgrid,mkqptequiv
!!      read_el_veloc,read_gkk
!!
!! CHILDREN
!!
!! SOURCE

integer function hkptrank_get_1kpt(krank,kpt) result(irank)

!Arguments ------------------------------------
!scalars
 type(hkptrank_t), intent(in) :: krank
!arrays
 double precision,intent(in) :: kpt(3)

!Local variables-------------------------------
!scalars
 character(len=500) :: msg
!arrays
 double precision :: redkpt(3),shift(3)

! *************************************************************************

 call wrap2_zero_one(kpt,redkpt,shift)

 irank = int(real(krank%max_linear_density)*(redkpt(1)+half+tol8 +&
&            real(krank%max_linear_density)*(redkpt(2)+half+tol8 +&
&            real(krank%max_linear_density)*(redkpt(3)+half+tol8))))

 if (irank > krank%max_rank .or. irank < 1) then
   write(msg,'(a,i0,a,i0)') ' rank should be between 1 and ', krank%max_rank, ' but got', irank
   MSG_ERROR(msg)
 end if

end function hkptrank_get_1kpt
!!***

!----------------------------------------------------------------------

!!****f* m_hkptrank/hkptrank_get_index
!!
!! NAME
!! hkptrank_get_index
!!
!! FUNCTION
!! This routine calculates the rank for one kpt
!!
!! INPUTS
!!  kpt = coordinates of kpoints
!!  krank = rank object for the k-grid we are using
!!
!! OUTPUT
!!  rank = rank of the kpoint
!!
!! PARENTS
!!      elphon,get_full_kgrid,integrate_gamma,integrate_gamma_alt,k_neighbors
!!      m_ddk,m_kptrank,m_nesting,m_pptools,m_tetrahedron,mkfskgrid,mkqptequiv
!!      read_el_veloc,read_gkk
!!
!! CHILDREN
!!
!! SOURCE

integer function hkptrank_get_index(krank,kpt) result(ikpt)

!Arguments ------------------------------------
!scalars
 type(hkptrank_t), intent(in) :: krank
!arrays
 double precision,intent(in) :: kpt(3)

!Local variables-------------------------------
!scalars
 integer :: ierr, irank

! *************************************************************************

 irank = hkptrank_get_1kpt(krank,kpt)
 call hashtable_get(krank%invrank,irank,ikpt,ierr)
 if (ierr/=0) MSG_ERROR('Did find mapping for k-point')

end function hkptrank_get_index
!!***

!----------------------------------------------------------------------

!!****f* m_hkptrank/hkptrank_free
!!
!! NAME
!! hkptrank_free
!!
!! FUNCTION
!! This routine deallocates the arrays in a kptrank_type structure
!!
!! INPUTS
!!  krank = object containing ranking and inverse ranking, to be deallocated
!!
!! PARENTS
!!      defs_elphon,get_full_kgrid,m_ddk,m_ebands,m_fstab,m_nesting,m_phgamma
!!      m_pptools,m_tetrahedron,mkfskgrid,mkqptequiv,order_fs_kpts,outelph
!!      read_el_veloc
!!
!! CHILDREN
!!
!! SOURCE

subroutine hkptrank_free (krank)

!Arguments ------------------------------------
!scalars
 type(hkptrank_t), intent(inout) :: krank
 call hashtable_free(krank%invrank)

end subroutine hkptrank_free
!!***

end module m_hkptrank
!!***

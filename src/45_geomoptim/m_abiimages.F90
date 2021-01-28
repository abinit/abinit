!!****m* ABINIT/m_abiimages
!! NAME
!! m_abiimages
!!
!! FUNCTION
!! This module contains definition one type to store the atomic
!! configuration of a set of images, and some methods to
!! create, manipulate and destroy the defined type
!!
!!
!! Datatypes :
!!
!! * abiimages     : Atomic configuration of a set of images
!!
!!
!!
!! Subroutines :
!!
!! * abiimages_ini
!! * abiimages_fin
!! * abiimages_bcast
!! * abiimages_compare
!!
!!
!! COPYRIGHT
!! Copyright (C) 2001-2020 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_abiimages

 use m_abicore
 use defs_basis
 use m_abihist
 use m_xmpi

 implicit none

 private

 public ::  abiimages
 public ::  abiimages_ini
 public ::  abiimages_fin


!----------------------------------------------------------------------
!!***

!!****t* m_abiimages/abiimages
!! NAME
!! abiimages
!!
!! FUNCTION
!! The type abiimages contains the configurations of N images,
!! The history of those images are store using N abihist objects
!! And the present configuration use a single abihist object
!!
!! It contains:
!!  img_present: An abihist object with the present values
!!               of the N images
!!  img_past:    An N-array of abihist objects containing the
!!               old values of the N images
!!
!! SOURCE

 type abiimages

! scalars
! Index of the last element on all records
    integer :: irecord
! Maximun size of the historical records
    integer :: nrecord
! Number of images
    integer :: nimage
! Number of images
    integer :: natom

! Atomic configurations of N images
    type(abihist)         :: img_present
! N historic records for N images
    type(abihist),allocatable :: img_past(:)


 end type abiimages
!!***

contains
!!***

!!****f* m_abiimages/abiimages_ini
!! NAME
!! abiimages_ini
!!
!! FUNCTION
!! Initialize the abiimages type
!!
!! INPUTS
!!  natom = Number of atoms per unitary cell for each image
!!  mxhist = Maximal number of records store per image
!!
!! OUTPUT
!!  abiimages <type(abiimages)> = The abiimages to initialize
!!
!! PARENTS
!!
!! CHILDREN
!!      abihist_fin
!!
!! SOURCE

subroutine abiimages_ini(images,nimages,natom,nrecord)

 implicit none

!Arguments ------------------------------------
 type(abiimages) ,intent(out) :: images
 integer         ,intent(in)  :: natom,nimages
 integer         ,intent(in)  :: nrecord
!Local variables-------------------------------
 integer :: i

! ***************************************************************

!Initialize indexes
 images%irecord=1
 images%nrecord=nrecord

!Allocate the present
 call abihist_ini(images%img_present,natom,nimages)

!Allocate the past
 ABI_MALLOC(images%img_past,(nimages))

 do i=1,nimages
    call abihist_ini(images%img_past(i),natom,nrecord)
 end do

end subroutine abiimages_ini
!!***

!----------------------------------------------------------------------

!!****f* m_abiimages/abiimages_fin
!! NAME
!! abiimages_fin
!!
!! FUNCTION
!! Deallocate all the pointers in a abiimages
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  images <type(abiimages)> = The abiimages to deallocate
!!
!! PARENTS
!!
!! CHILDREN
!!      abihist_fin
!!
!! NOTES
!!
!! SOURCE

subroutine abiimages_fin(images)

 implicit none

!Arguments ------------------------------------
 type(abiimages),intent(inout) :: images
!Local variables-------------------------------
 integer :: i

! ***************************************************************

 call abihist_fin(images%img_present)

!Nullify the past
 do i=1,images%natom
    call abihist_fin(images%img_past(i))
 end do

 ABI_FREE(images%img_past)

end subroutine abiimages_fin
!!***

end module m_abiimages
!!***

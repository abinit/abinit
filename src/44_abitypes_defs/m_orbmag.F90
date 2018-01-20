!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_orbmag
!! NAME
!!  m_orbmag
!!
!! FUNCTION
!!  This module contains the declaration of data types and methods
!!  used to handle orbital magnetization
!!
!! COPYRIGHT
!! Copyright (C) 2011-2017 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
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

module m_orbmag

 use defs_basis
 use m_profiling_abi
 use m_errors

 use m_pawcprj, only : pawcprj_type, pawcprj_free

 implicit none

 private
!!***


!!****t* defs_datatypes/orbmag_type
!! NAME
!! orbmag_type
!!
!! FUNCTION
!! variables used in orbital magnetism calculation
!!
!! SOURCE

 type, public :: orbmag_type

! WARNING : if you modify this datatype, please check whether there might be creation/destruction/copy routines,
! declared in another part of ABINIT, that might need to take into account your modification.

! Integer variables
  integer :: orbmag              ! value of orbmag input variable in use
  integer :: fmkmem              ! number of k-points in the FBZ per cpu
  integer :: fmkmem_max          ! max of fmkmem
  integer :: fnkpt               ! number of k-points in the FBZ
  integer :: lmax
  integer :: lmnmax
  integer :: lmn2max
  integer :: mkmem_max           ! max of mkmem
  integer :: natom               ! number of atoms in unit cell
  integer :: my_natom            ! number of atoms treated by current proc
  integer :: mband_occ           ! max number of occupied bands (over spin)
                                 ! this number must be the same for every k
  integer :: nspinor             ! nspinor input from data set
  integer :: nsym
  integer :: usepaw              ! 1 if a PAW calculation, 0 else

! Real(dp) scalars
  real(dp) :: sdeg               ! spin degeneracy: sdeg = 2 if nsppol = 1

  ! Real(dp) arrays
  real(dp) :: chern(2,3)           ! result of chern number calculation
  
  real(dp) :: dkvecs(3,3)        ! dkvec(:,idir) = vector between a k-poinit
                                 ! and its nearest neighbour along idir
! Integer pointers
  integer, allocatable :: atom_indsym(:,:,:) ! atom_indsym(4,nsym,natom)
                                         ! this is data on how the symmetries map the atoms in the cell
                                         ! see symatm.F90 for full description
  integer, allocatable :: cgindex(:,:)    ! cgindex(nkpt,nsppol)
                                      ! for each k-point, stores the location
                                      ! of the WF in the cg array
  integer, allocatable :: cprjindex(:,:)  ! cprjindex(nkpt,nsppol)
                                      ! for each k-point, stores the location
                                      ! of the cprj in the cprj array (used only
                                      ! for PAW calculations)
  integer, allocatable :: fkgindex(:)     ! same as kgindex, but defined
                                      ! for the FBZ and intended to use
                                      ! with pwindf
  integer, allocatable :: ikpt_dk(:,:,:)  ! ikpt_dk(nkpt,2,3)
                                      ! ikpt_dp(ikpt,ii,idir) = index of the
                                      ! k-point at k+dk (ii=1) and k-dk (ii=2)
  integer, allocatable :: indkk_f2ibz(:,:)   ! indkk_f2ibz(1:dtorbmag%fnkpt,1:6)
                                         ! information needed to fold a
                                         ! k-point in the FBZ into the IBZ;
                                         ! the second index (1:6)
                                         ! is as described in listkk
  integer, allocatable :: i2fbz(:)           ! i2fbz(1:nkpt) gives index of IBZ
                                         ! k-points in the FBZ k-point list

  integer, allocatable :: kgindex(:)      ! kgind(nkpt)
                                      ! kgind(ikpt) = ikg

  integer, allocatable :: lmn_size(:)        ! lmn_size(ntypat)
  integer, allocatable :: lmn2_size(:)       ! lmn2_size(ntypat)

  integer, allocatable :: nband_occ(:)       ! nband_occ(nsppol) = actual number of occupied bands
                                             !  can be different for spin up and down!!!
! Real(dp) allocatables

  real(dp), allocatable :: fkptns(:,:)       ! fkptns(3,1:dtorbmag%fnkpt)
                                         ! k-points in FBZ

  real(dp), allocatable :: zarot(:,:,:,:)
   !  zarot(l_size_max,l_size_max,l_max,nsym)
   !  Coeffs of the transformation of real spherical
   !  harmonics under the symmetry operations. These are needed when the
   ! cprj's need to be computed in the full BZ, that is,
   ! in the PAW case with kptopt /= 3.

 end type orbmag_type

 ! Bound methods:
 public :: destroy_orbmag
!!***

contains

!!****f* m_orbmag/destroy_orbmag
!! NAME
!!
!! FUNCTION
!!   deallocate fields in orbmag structure
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine destroy_orbmag(dtorbmag)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'destroy_orbmag'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!array
 type(orbmag_type),intent(inout) :: dtorbmag

! ************************************************************************

! Integer pointers
  if(allocated(dtorbmag%atom_indsym))  then
    ABI_DEALLOCATE(dtorbmag%atom_indsym)
  end if
  if(allocated(dtorbmag%cgindex))  then
    ABI_DEALLOCATE(dtorbmag%cgindex)
  end if
  if(allocated(dtorbmag%cprjindex))  then
    ABI_DEALLOCATE(dtorbmag%cprjindex)
  end if
  if(allocated(dtorbmag%fkgindex))  then
    ABI_DEALLOCATE(dtorbmag%fkgindex)
  end if
  if(allocated(dtorbmag%ikpt_dk))  then
    ABI_DEALLOCATE(dtorbmag%ikpt_dk)
  end if
  if(allocated(dtorbmag%indkk_f2ibz))  then
    ABI_DEALLOCATE(dtorbmag%indkk_f2ibz)
  end if
  if(allocated(dtorbmag%i2fbz))  then
    ABI_DEALLOCATE(dtorbmag%i2fbz)
  end if
  if(allocated(dtorbmag%kgindex))  then
    ABI_DEALLOCATE(dtorbmag%kgindex)
  end if
  if(allocated(dtorbmag%lmn_size))  then
    ABI_DEALLOCATE(dtorbmag%lmn_size)
  end if
  if(allocated(dtorbmag%lmn2_size))  then
    ABI_DEALLOCATE(dtorbmag%lmn2_size)
  end if
  if(allocated(dtorbmag%nband_occ))  then
    ABI_DEALLOCATE(dtorbmag%nband_occ)
  end if
! Real(dp) pointers

  if(allocated(dtorbmag%fkptns))  then
    ABI_DEALLOCATE(dtorbmag%fkptns)
  end if
  if(allocated(dtorbmag%zarot))  then
    ABI_DEALLOCATE(dtorbmag%zarot)
  end if

end subroutine destroy_orbmag
!!***

end module m_orbmag
!!***

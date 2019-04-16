!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_htetrahedron
!! NAME
!! m_htetrahedron
!!
!! FUNCTION
!!  module for tetrahedron interpolation of DOS and similar quantities
!!  depends on m_kpt_rank.
!!  Uses some functions from a previous implementation by MJV
!!  The new implementation if based on spglib and kpclib by Atsushi Togo
!!  after a discussion with in on the APS 2019 where he provided
!!  details of his implementation.
!!
!! COPYRIGHT
!!  Copyright (C) 2010-2019 ABINIT group (HM,MJV)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! TODO
!!  1) Test carefully the case of degenerate tethraedron
!!  2) Add options to get only delta and/or theta ?
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

module m_htetrahedron

 use defs_basis
 use m_abicore
 use m_kptrank
 use m_simtet,          only : sim0onei, SIM0TWOI
 use m_xmpi

implicit none

private
!!***

integer, parameter :: TETRA_SIZE = 24
integer, parameter :: TETRA_STEP = 24

!!****t* m_htetrahedron/t_htetra_bucket
!! NAME
!! t_htetra_bucket
!!
!! FUNCTION
!! Store a bunch of tetrahedra
!!
!! SOURCE

type :: htetra_bucket

  integer,pointer :: indexes(:,:)

end type htetra_bucket

!!****t* m_htetrahedron/t_htetrap
!! NAME
!! t_htetrap
!!
!! FUNCTION
!! Pointer to tetrahedra
!!
!! SOURCE

type :: tetrap

  integer,pointer :: p(:)

end type tetrap

!!****t* m_htetrahedron/t_htetrak
!! NAME
!! t_htetrak
!!
!! FUNCTION
!! Pointer to tetrahedra associated to a k-point
!!
!! SOURCE

type :: htetrak

  integer :: tetra_count
  integer :: tetra_total
  type(tetrap),allocatable :: tetra(:)

end type htetrak


!!****t* m_htetrahedron/t_htetrahedron
!! NAME
!! t_htetrahedron
!!
!! FUNCTION
!! tetrahedron geometry object
!!
!! SOURCE

type, public :: t_htetrahedron

  integer :: nkibz
  ! Number of points in the irreducible Brillouin zone

  integer :: nkbz
  ! Number of points in the full Brillouin zone

  integer :: nunique_tetra
  ! Number of unique tetrahedron

  integer :: nibz_tetra
  ! Number of ibz tetrahedron

  real(dp)  :: vv
  ! volume of the tetrahedra

  real(dp) :: klatt(3, 3)
  ! reciprocal of lattice vectors for full kpoint grid

  real(dp),allocatable :: ibz_weights(:)
  ! weights for integration

  type(htetrak),allocatable :: ibz(:)
  ! indexes of the tetrahedra for each k-point

  type(htetra_bucket),allocatable :: unique_tetra(:)
  ! indexes of the unique tetrahedra

end type t_htetrahedron

public :: htetra_init            ! Initialize the object
public :: htetra_free            ! Free memory
public :: htetra_print           ! Print information about tetrahedron object
public :: htetra_get_onewk       ! Calculate integration weights and their derivatives for a single k-point in the IBZ.
public :: htetra_get_onewk_wvals ! Similar to tetra_get_onewk_wvals but receives arbitrary list of frequency points.
public :: htetra_get_onewk_wvals_zinv ! Calculate integration weights for 1/(z-E(k)) for a single k-point in the IBZ.
public :: htetra_blochl_weights  ! And interface to help to facilitate the transition to the new tetrahedron implementation
!!***

contains
!!***

!----------------------------------------------------------------------

!!****f* m_htetrahedron/htetra_init
!! NAME
!! htetra_init
!!
!! FUNCTION
!! get tetrahedra characterized by apexes
!!
!! INPUTS
!!  bz2ibz(nkpt_fullbz)=indexes of irred kpoints equivalent to kpt_fullbz
!!  gprimd(3,3) = reciprocal space vectors
!!  klatt(3,3)=reciprocal of lattice vectors for full kpoint grid
!!  kpt_fullbz(3,nkpt_fullbz)=kpoints in full brillouin zone
!!  nkpt_fullbz=number of kpoints in full brillouin zone
!!  comm= MPI communicator
!!
!! OUTPUT
!!  tetra%ibz(4,24,nkibz)=for each k-point, the indexes in the IBZ
!!  tetra%vv = tetrahedron volume divided by full BZ volume
!!
!! PARENTS
!!      ep_el_weights,ep_fs_weights,ep_ph_weights,m_fstab,m_kpts,m_phonons
!!      thmeig
!!
!! CHILDREN
!!
!! SOURCE

subroutine htetra_init(tetra, bz2ibz, gprimd, klatt, kpt_fullbz, nkpt_fullbz, kpt_ibz, nkpt_ibz, ierr, errorstring, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt_fullbz, nkpt_ibz, comm
 integer, intent(out) :: ierr
 character(len=80), intent(out) :: errorstring
 type(t_htetrahedron),intent(out) :: tetra
!arrays
 integer,intent(in) :: bz2ibz(nkpt_fullbz)
 real(dp) ,intent(in) :: gprimd(3,3),klatt(3,3),kpt_fullbz(3,nkpt_fullbz),kpt_ibz(3,nkpt_ibz)

!Local variables-------------------------------
!scalars
 type(kptrank_type) :: kptrank_t
 integer :: ikpt2,isummit,itetra,jtetra, max_tetra_count
 integer :: ikibz,ikbz,idiag,ihash,min_idiag,my_rank,nprocs
 integer :: symrankkpt, max_ntetra, tetra_total, total_ntetra, ntetra, hash
 real(dp) :: rcvol,length,min_length
 character(len=500) :: msg
!arrays
 integer,pointer :: indexes(:,:)
 integer :: tetra_ibz(4), tetra_mibz(0:4), tetra_count(nkpt_ibz)
 integer :: tetra_shifts(3,4,24,4)  ! 3 dimensions, 4 summits, 24 tetrahedra, 4 main diagonals
 integer :: tetra_shifts_6(3,4,6,1) ! 3 dimensions, 4 summits, 6 tetrahedra, 4 main diagonals
 integer :: main_diagonals(3,4)
 real(dp)  :: k1(3),k2(3),k3(3),diag(3)

! *********************************************************************

 ! Use the shifts from kpclib developed by Atsushi Togo
 ! This part is produced by a python script
 tetra_shifts(:, 1, 1,1) = [  0,  0,  0]
 tetra_shifts(:, 2, 1,1) = [  1,  0,  0]
 tetra_shifts(:, 3, 1,1) = [  1,  1,  0]
 tetra_shifts(:, 4, 1,1) = [  1,  1,  1]
 tetra_shifts(:, 1, 2,1) = [  0,  0,  0]
 tetra_shifts(:, 2, 2,1) = [  1,  0,  0]
 tetra_shifts(:, 3, 2,1) = [  1,  0,  1]
 tetra_shifts(:, 4, 2,1) = [  1,  1,  1]
 tetra_shifts(:, 1, 3,1) = [  0,  0,  0]
 tetra_shifts(:, 2, 3,1) = [  0,  1,  0]
 tetra_shifts(:, 3, 3,1) = [  1,  1,  0]
 tetra_shifts(:, 4, 3,1) = [  1,  1,  1]
 tetra_shifts(:, 1, 4,1) = [  0,  0,  0]
 tetra_shifts(:, 2, 4,1) = [  0,  1,  0]
 tetra_shifts(:, 3, 4,1) = [  0,  1,  1]
 tetra_shifts(:, 4, 4,1) = [  1,  1,  1]
 tetra_shifts(:, 1, 5,1) = [  0,  0,  0]
 tetra_shifts(:, 2, 5,1) = [  0,  0,  1]
 tetra_shifts(:, 3, 5,1) = [  1,  0,  1]
 tetra_shifts(:, 4, 5,1) = [  1,  1,  1]
 tetra_shifts(:, 1, 6,1) = [  0,  0,  0]
 tetra_shifts(:, 2, 6,1) = [  0,  0,  1]
 tetra_shifts(:, 3, 6,1) = [  0,  1,  1]
 tetra_shifts(:, 4, 6,1) = [  1,  1,  1]
 tetra_shifts(:, 1, 7,1) = [  0,  0,  0]
 tetra_shifts(:, 2, 7,1) = [  0,  1,  0]
 tetra_shifts(:, 3, 7,1) = [  0,  1,  1]
 tetra_shifts(:, 4, 7,1) = [ -1,  0,  0]
 tetra_shifts(:, 1, 8,1) = [  0,  0,  0]
 tetra_shifts(:, 2, 8,1) = [  0,  0,  1]
 tetra_shifts(:, 3, 8,1) = [  0,  1,  1]
 tetra_shifts(:, 4, 8,1) = [ -1,  0,  0]
 tetra_shifts(:, 1, 9,1) = [  0,  0,  0]
 tetra_shifts(:, 2, 9,1) = [  1,  0,  0]
 tetra_shifts(:, 3, 9,1) = [  1,  0,  1]
 tetra_shifts(:, 4, 9,1) = [  0, -1,  0]
 tetra_shifts(:, 1,10,1) = [  0,  0,  0]
 tetra_shifts(:, 2,10,1) = [  0,  0,  1]
 tetra_shifts(:, 3,10,1) = [  1,  0,  1]
 tetra_shifts(:, 4,10,1) = [  0, -1,  0]
 tetra_shifts(:, 1,11,1) = [  0,  0,  0]
 tetra_shifts(:, 2,11,1) = [  0,  0,  1]
 tetra_shifts(:, 3,11,1) = [ -1, -1,  0]
 tetra_shifts(:, 4,11,1) = [  0, -1,  0]
 tetra_shifts(:, 1,12,1) = [  0,  0,  0]
 tetra_shifts(:, 2,12,1) = [  0,  0,  1]
 tetra_shifts(:, 3,12,1) = [ -1, -1,  0]
 tetra_shifts(:, 4,12,1) = [ -1,  0,  0]
 tetra_shifts(:, 1,13,1) = [  0,  0,  0]
 tetra_shifts(:, 2,13,1) = [  1,  0,  0]
 tetra_shifts(:, 3,13,1) = [  1,  1,  0]
 tetra_shifts(:, 4,13,1) = [  0,  0, -1]
 tetra_shifts(:, 1,14,1) = [  0,  0,  0]
 tetra_shifts(:, 2,14,1) = [  0,  1,  0]
 tetra_shifts(:, 3,14,1) = [  1,  1,  0]
 tetra_shifts(:, 4,14,1) = [  0,  0, -1]
 tetra_shifts(:, 1,15,1) = [  0,  0,  0]
 tetra_shifts(:, 2,15,1) = [  0,  1,  0]
 tetra_shifts(:, 3,15,1) = [ -1,  0, -1]
 tetra_shifts(:, 4,15,1) = [  0,  0, -1]
 tetra_shifts(:, 1,16,1) = [  0,  0,  0]
 tetra_shifts(:, 2,16,1) = [  0,  1,  0]
 tetra_shifts(:, 3,16,1) = [ -1,  0, -1]
 tetra_shifts(:, 4,16,1) = [ -1,  0,  0]
 tetra_shifts(:, 1,17,1) = [  0,  0,  0]
 tetra_shifts(:, 2,17,1) = [  1,  0,  0]
 tetra_shifts(:, 3,17,1) = [  0, -1, -1]
 tetra_shifts(:, 4,17,1) = [  0,  0, -1]
 tetra_shifts(:, 1,18,1) = [  0,  0,  0]
 tetra_shifts(:, 2,18,1) = [  1,  0,  0]
 tetra_shifts(:, 3,18,1) = [  0, -1, -1]
 tetra_shifts(:, 4,18,1) = [  0, -1,  0]
 tetra_shifts(:, 1,19,1) = [  0,  0,  0]
 tetra_shifts(:, 2,19,1) = [ -1, -1, -1]
 tetra_shifts(:, 3,19,1) = [  0, -1, -1]
 tetra_shifts(:, 4,19,1) = [  0,  0, -1]
 tetra_shifts(:, 1,20,1) = [  0,  0,  0]
 tetra_shifts(:, 2,20,1) = [ -1, -1, -1]
 tetra_shifts(:, 3,20,1) = [  0, -1, -1]
 tetra_shifts(:, 4,20,1) = [  0, -1,  0]
 tetra_shifts(:, 1,21,1) = [  0,  0,  0]
 tetra_shifts(:, 2,21,1) = [ -1, -1, -1]
 tetra_shifts(:, 3,21,1) = [ -1,  0, -1]
 tetra_shifts(:, 4,21,1) = [  0,  0, -1]
 tetra_shifts(:, 1,22,1) = [  0,  0,  0]
 tetra_shifts(:, 2,22,1) = [ -1, -1, -1]
 tetra_shifts(:, 3,22,1) = [ -1,  0, -1]
 tetra_shifts(:, 4,22,1) = [ -1,  0,  0]
 tetra_shifts(:, 1,23,1) = [  0,  0,  0]
 tetra_shifts(:, 2,23,1) = [ -1, -1, -1]
 tetra_shifts(:, 3,23,1) = [ -1, -1,  0]
 tetra_shifts(:, 4,23,1) = [  0, -1,  0]
 tetra_shifts(:, 1,24,1) = [  0,  0,  0]
 tetra_shifts(:, 2,24,1) = [ -1, -1, -1]
 tetra_shifts(:, 3,24,1) = [ -1, -1,  0]
 tetra_shifts(:, 4,24,1) = [ -1,  0,  0]
 tetra_shifts(:, 1, 1,2) = [  0,  0,  0]
 tetra_shifts(:, 2, 1,2) = [  1,  0,  0]
 tetra_shifts(:, 3, 1,2) = [  0,  1,  0]
 tetra_shifts(:, 4, 1,2) = [  0,  1,  1]
 tetra_shifts(:, 1, 2,2) = [  0,  0,  0]
 tetra_shifts(:, 2, 2,2) = [  1,  0,  0]
 tetra_shifts(:, 3, 2,2) = [  0,  0,  1]
 tetra_shifts(:, 4, 2,2) = [  0,  1,  1]
 tetra_shifts(:, 1, 3,2) = [  0,  0,  0]
 tetra_shifts(:, 2, 3,2) = [ -1,  1,  0]
 tetra_shifts(:, 3, 3,2) = [ -1,  1,  1]
 tetra_shifts(:, 4, 3,2) = [ -1,  0,  0]
 tetra_shifts(:, 1, 4,2) = [  0,  0,  0]
 tetra_shifts(:, 2, 4,2) = [ -1,  0,  1]
 tetra_shifts(:, 3, 4,2) = [ -1,  1,  1]
 tetra_shifts(:, 4, 4,2) = [ -1,  0,  0]
 tetra_shifts(:, 1, 5,2) = [  0,  0,  0]
 tetra_shifts(:, 2, 5,2) = [ -1,  1,  0]
 tetra_shifts(:, 3, 5,2) = [  0,  1,  0]
 tetra_shifts(:, 4, 5,2) = [ -1,  1,  1]
 tetra_shifts(:, 1, 6,2) = [  0,  0,  0]
 tetra_shifts(:, 2, 6,2) = [  0,  1,  0]
 tetra_shifts(:, 3, 6,2) = [ -1,  1,  1]
 tetra_shifts(:, 4, 6,2) = [  0,  1,  1]
 tetra_shifts(:, 1, 7,2) = [  0,  0,  0]
 tetra_shifts(:, 2, 7,2) = [ -1,  0,  1]
 tetra_shifts(:, 3, 7,2) = [  0,  0,  1]
 tetra_shifts(:, 4, 7,2) = [ -1,  1,  1]
 tetra_shifts(:, 1, 8,2) = [  0,  0,  0]
 tetra_shifts(:, 2, 8,2) = [  0,  0,  1]
 tetra_shifts(:, 3, 8,2) = [ -1,  1,  1]
 tetra_shifts(:, 4, 8,2) = [  0,  1,  1]
 tetra_shifts(:, 1, 9,2) = [  0,  0,  0]
 tetra_shifts(:, 2, 9,2) = [  0,  0,  1]
 tetra_shifts(:, 3, 9,2) = [  0, -1,  0]
 tetra_shifts(:, 4, 9,2) = [  1, -1,  0]
 tetra_shifts(:, 1,10,2) = [  0,  0,  0]
 tetra_shifts(:, 2,10,2) = [  1,  0,  0]
 tetra_shifts(:, 3,10,2) = [  0,  0,  1]
 tetra_shifts(:, 4,10,2) = [  1, -1,  0]
 tetra_shifts(:, 1,11,2) = [  0,  0,  0]
 tetra_shifts(:, 2,11,2) = [ -1,  0,  1]
 tetra_shifts(:, 3,11,2) = [  0, -1,  0]
 tetra_shifts(:, 4,11,2) = [ -1,  0,  0]
 tetra_shifts(:, 1,12,2) = [  0,  0,  0]
 tetra_shifts(:, 2,12,2) = [ -1,  0,  1]
 tetra_shifts(:, 3,12,2) = [  0,  0,  1]
 tetra_shifts(:, 4,12,2) = [  0, -1,  0]
 tetra_shifts(:, 1,13,2) = [  0,  0,  0]
 tetra_shifts(:, 2,13,2) = [  0,  1,  0]
 tetra_shifts(:, 3,13,2) = [  0,  0, -1]
 tetra_shifts(:, 4,13,2) = [  1,  0, -1]
 tetra_shifts(:, 1,14,2) = [  0,  0,  0]
 tetra_shifts(:, 2,14,2) = [  1,  0,  0]
 tetra_shifts(:, 3,14,2) = [  0,  1,  0]
 tetra_shifts(:, 4,14,2) = [  1,  0, -1]
 tetra_shifts(:, 1,15,2) = [  0,  0,  0]
 tetra_shifts(:, 2,15,2) = [ -1,  1,  0]
 tetra_shifts(:, 3,15,2) = [  0,  0, -1]
 tetra_shifts(:, 4,15,2) = [ -1,  0,  0]
 tetra_shifts(:, 1,16,2) = [  0,  0,  0]
 tetra_shifts(:, 2,16,2) = [ -1,  1,  0]
 tetra_shifts(:, 3,16,2) = [  0,  1,  0]
 tetra_shifts(:, 4,16,2) = [  0,  0, -1]
 tetra_shifts(:, 1,17,2) = [  0,  0,  0]
 tetra_shifts(:, 2,17,2) = [  0, -1, -1]
 tetra_shifts(:, 3,17,2) = [  1, -1, -1]
 tetra_shifts(:, 4,17,2) = [  0,  0, -1]
 tetra_shifts(:, 1,18,2) = [  0,  0,  0]
 tetra_shifts(:, 2,18,2) = [  0, -1, -1]
 tetra_shifts(:, 3,18,2) = [  1, -1, -1]
 tetra_shifts(:, 4,18,2) = [  0, -1,  0]
 tetra_shifts(:, 1,19,2) = [  0,  0,  0]
 tetra_shifts(:, 2,19,2) = [  1, -1, -1]
 tetra_shifts(:, 3,19,2) = [  0,  0, -1]
 tetra_shifts(:, 4,19,2) = [  1,  0, -1]
 tetra_shifts(:, 1,20,2) = [  0,  0,  0]
 tetra_shifts(:, 2,20,2) = [  1,  0,  0]
 tetra_shifts(:, 3,20,2) = [  1, -1, -1]
 tetra_shifts(:, 4,20,2) = [  1,  0, -1]
 tetra_shifts(:, 1,21,2) = [  0,  0,  0]
 tetra_shifts(:, 2,21,2) = [  1, -1, -1]
 tetra_shifts(:, 3,21,2) = [  0, -1,  0]
 tetra_shifts(:, 4,21,2) = [  1, -1,  0]
 tetra_shifts(:, 1,22,2) = [  0,  0,  0]
 tetra_shifts(:, 2,22,2) = [  1,  0,  0]
 tetra_shifts(:, 3,22,2) = [  1, -1, -1]
 tetra_shifts(:, 4,22,2) = [  1, -1,  0]
 tetra_shifts(:, 1,23,2) = [  0,  0,  0]
 tetra_shifts(:, 2,23,2) = [  0, -1, -1]
 tetra_shifts(:, 3,23,2) = [  0,  0, -1]
 tetra_shifts(:, 4,23,2) = [ -1,  0,  0]
 tetra_shifts(:, 1,24,2) = [  0,  0,  0]
 tetra_shifts(:, 2,24,2) = [  0, -1, -1]
 tetra_shifts(:, 3,24,2) = [  0, -1,  0]
 tetra_shifts(:, 4,24,2) = [ -1,  0,  0]
 tetra_shifts(:, 1, 1,3) = [  0,  0,  0]
 tetra_shifts(:, 2, 1,3) = [  1,  0,  0]
 tetra_shifts(:, 3, 1,3) = [  0,  1,  0]
 tetra_shifts(:, 4, 1,3) = [  1,  0,  1]
 tetra_shifts(:, 1, 2,3) = [  0,  0,  0]
 tetra_shifts(:, 2, 2,3) = [  0,  1,  0]
 tetra_shifts(:, 3, 2,3) = [  0,  0,  1]
 tetra_shifts(:, 4, 2,3) = [  1,  0,  1]
 tetra_shifts(:, 1, 3,3) = [  0,  0,  0]
 tetra_shifts(:, 2, 3,3) = [ -1,  1,  0]
 tetra_shifts(:, 3, 3,3) = [  0,  0,  1]
 tetra_shifts(:, 4, 3,3) = [ -1,  0,  0]
 tetra_shifts(:, 1, 4,3) = [  0,  0,  0]
 tetra_shifts(:, 2, 4,3) = [ -1,  1,  0]
 tetra_shifts(:, 3, 4,3) = [  0,  1,  0]
 tetra_shifts(:, 4, 4,3) = [  0,  0,  1]
 tetra_shifts(:, 1, 5,3) = [  0,  0,  0]
 tetra_shifts(:, 2, 5,3) = [  1, -1,  1]
 tetra_shifts(:, 3, 5,3) = [  0, -1,  0]
 tetra_shifts(:, 4, 5,3) = [  1, -1,  0]
 tetra_shifts(:, 1, 6,3) = [  0,  0,  0]
 tetra_shifts(:, 2, 6,3) = [  0, -1,  1]
 tetra_shifts(:, 3, 6,3) = [  1, -1,  1]
 tetra_shifts(:, 4, 6,3) = [  0, -1,  0]
 tetra_shifts(:, 1, 7,3) = [  0,  0,  0]
 tetra_shifts(:, 2, 7,3) = [  1,  0,  0]
 tetra_shifts(:, 3, 7,3) = [  1, -1,  1]
 tetra_shifts(:, 4, 7,3) = [  1, -1,  0]
 tetra_shifts(:, 1, 8,3) = [  0,  0,  0]
 tetra_shifts(:, 2, 8,3) = [  1,  0,  0]
 tetra_shifts(:, 3, 8,3) = [  1, -1,  1]
 tetra_shifts(:, 4, 8,3) = [  1,  0,  1]
 tetra_shifts(:, 1, 9,3) = [  0,  0,  0]
 tetra_shifts(:, 2, 9,3) = [  0, -1,  1]
 tetra_shifts(:, 3, 9,3) = [  1, -1,  1]
 tetra_shifts(:, 4, 9,3) = [  0,  0,  1]
 tetra_shifts(:, 1,10,3) = [  0,  0,  0]
 tetra_shifts(:, 2,10,3) = [  1, -1,  1]
 tetra_shifts(:, 3,10,3) = [  0,  0,  1]
 tetra_shifts(:, 4,10,3) = [  1,  0,  1]
 tetra_shifts(:, 1,11,3) = [  0,  0,  0]
 tetra_shifts(:, 2,11,3) = [  0, -1,  1]
 tetra_shifts(:, 3,11,3) = [  0, -1,  0]
 tetra_shifts(:, 4,11,3) = [ -1,  0,  0]
 tetra_shifts(:, 1,12,3) = [  0,  0,  0]
 tetra_shifts(:, 2,12,3) = [  0, -1,  1]
 tetra_shifts(:, 3,12,3) = [  0,  0,  1]
 tetra_shifts(:, 4,12,3) = [ -1,  0,  0]
 tetra_shifts(:, 1,13,3) = [  0,  0,  0]
 tetra_shifts(:, 2,13,3) = [  1,  0,  0]
 tetra_shifts(:, 3,13,3) = [  0,  0, -1]
 tetra_shifts(:, 4,13,3) = [  0,  1, -1]
 tetra_shifts(:, 1,14,3) = [  0,  0,  0]
 tetra_shifts(:, 2,14,3) = [  1,  0,  0]
 tetra_shifts(:, 3,14,3) = [  0,  1,  0]
 tetra_shifts(:, 4,14,3) = [  0,  1, -1]
 tetra_shifts(:, 1,15,3) = [  0,  0,  0]
 tetra_shifts(:, 2,15,3) = [ -1,  0, -1]
 tetra_shifts(:, 3,15,3) = [  0,  0, -1]
 tetra_shifts(:, 4,15,3) = [ -1,  1, -1]
 tetra_shifts(:, 1,16,3) = [  0,  0,  0]
 tetra_shifts(:, 2,16,3) = [ -1,  0, -1]
 tetra_shifts(:, 3,16,3) = [ -1,  1, -1]
 tetra_shifts(:, 4,16,3) = [ -1,  0,  0]
 tetra_shifts(:, 1,17,3) = [  0,  0,  0]
 tetra_shifts(:, 2,17,3) = [  0,  0, -1]
 tetra_shifts(:, 3,17,3) = [ -1,  1, -1]
 tetra_shifts(:, 4,17,3) = [  0,  1, -1]
 tetra_shifts(:, 1,18,3) = [  0,  0,  0]
 tetra_shifts(:, 2,18,3) = [  0,  1,  0]
 tetra_shifts(:, 3,18,3) = [ -1,  1, -1]
 tetra_shifts(:, 4,18,3) = [  0,  1, -1]
 tetra_shifts(:, 1,19,3) = [  0,  0,  0]
 tetra_shifts(:, 2,19,3) = [ -1,  1,  0]
 tetra_shifts(:, 3,19,3) = [ -1,  1, -1]
 tetra_shifts(:, 4,19,3) = [ -1,  0,  0]
 tetra_shifts(:, 1,20,3) = [  0,  0,  0]
 tetra_shifts(:, 2,20,3) = [ -1,  1,  0]
 tetra_shifts(:, 3,20,3) = [  0,  1,  0]
 tetra_shifts(:, 4,20,3) = [ -1,  1, -1]
 tetra_shifts(:, 1,21,3) = [  0,  0,  0]
 tetra_shifts(:, 2,21,3) = [  0,  0, -1]
 tetra_shifts(:, 3,21,3) = [  0, -1,  0]
 tetra_shifts(:, 4,21,3) = [  1, -1,  0]
 tetra_shifts(:, 1,22,3) = [  0,  0,  0]
 tetra_shifts(:, 2,22,3) = [  1,  0,  0]
 tetra_shifts(:, 3,22,3) = [  0,  0, -1]
 tetra_shifts(:, 4,22,3) = [  1, -1,  0]
 tetra_shifts(:, 1,23,3) = [  0,  0,  0]
 tetra_shifts(:, 2,23,3) = [ -1,  0, -1]
 tetra_shifts(:, 3,23,3) = [  0,  0, -1]
 tetra_shifts(:, 4,23,3) = [  0, -1,  0]
 tetra_shifts(:, 1,24,3) = [  0,  0,  0]
 tetra_shifts(:, 2,24,3) = [ -1,  0, -1]
 tetra_shifts(:, 3,24,3) = [  0, -1,  0]
 tetra_shifts(:, 4,24,3) = [ -1,  0,  0]
 tetra_shifts(:, 1, 1,4) = [  0,  0,  0]
 tetra_shifts(:, 2, 1,4) = [  1,  0,  0]
 tetra_shifts(:, 3, 1,4) = [  1,  1,  0]
 tetra_shifts(:, 4, 1,4) = [  0,  0,  1]
 tetra_shifts(:, 1, 2,4) = [  0,  0,  0]
 tetra_shifts(:, 2, 2,4) = [  0,  1,  0]
 tetra_shifts(:, 3, 2,4) = [  1,  1,  0]
 tetra_shifts(:, 4, 2,4) = [  0,  0,  1]
 tetra_shifts(:, 1, 3,4) = [  0,  0,  0]
 tetra_shifts(:, 2, 3,4) = [  0,  1,  0]
 tetra_shifts(:, 3, 3,4) = [ -1,  0,  1]
 tetra_shifts(:, 4, 3,4) = [ -1,  0,  0]
 tetra_shifts(:, 1, 4,4) = [  0,  0,  0]
 tetra_shifts(:, 2, 4,4) = [  0,  1,  0]
 tetra_shifts(:, 3, 4,4) = [ -1,  0,  1]
 tetra_shifts(:, 4, 4,4) = [  0,  0,  1]
 tetra_shifts(:, 1, 5,4) = [  0,  0,  0]
 tetra_shifts(:, 2, 5,4) = [  1,  0,  0]
 tetra_shifts(:, 3, 5,4) = [  0, -1,  1]
 tetra_shifts(:, 4, 5,4) = [  0, -1,  0]
 tetra_shifts(:, 1, 6,4) = [  0,  0,  0]
 tetra_shifts(:, 2, 6,4) = [  1,  0,  0]
 tetra_shifts(:, 3, 6,4) = [  0, -1,  1]
 tetra_shifts(:, 4, 6,4) = [  0,  0,  1]
 tetra_shifts(:, 1, 7,4) = [  0,  0,  0]
 tetra_shifts(:, 2, 7,4) = [ -1, -1,  1]
 tetra_shifts(:, 3, 7,4) = [ -1, -1,  0]
 tetra_shifts(:, 4, 7,4) = [  0, -1,  0]
 tetra_shifts(:, 1, 8,4) = [  0,  0,  0]
 tetra_shifts(:, 2, 8,4) = [ -1, -1,  1]
 tetra_shifts(:, 3, 8,4) = [ -1, -1,  0]
 tetra_shifts(:, 4, 8,4) = [ -1,  0,  0]
 tetra_shifts(:, 1, 9,4) = [  0,  0,  0]
 tetra_shifts(:, 2, 9,4) = [ -1, -1,  1]
 tetra_shifts(:, 3, 9,4) = [  0, -1,  1]
 tetra_shifts(:, 4, 9,4) = [  0, -1,  0]
 tetra_shifts(:, 1,10,4) = [  0,  0,  0]
 tetra_shifts(:, 2,10,4) = [ -1, -1,  1]
 tetra_shifts(:, 3,10,4) = [ -1,  0,  1]
 tetra_shifts(:, 4,10,4) = [ -1,  0,  0]
 tetra_shifts(:, 1,11,4) = [  0,  0,  0]
 tetra_shifts(:, 2,11,4) = [ -1, -1,  1]
 tetra_shifts(:, 3,11,4) = [  0, -1,  1]
 tetra_shifts(:, 4,11,4) = [  0,  0,  1]
 tetra_shifts(:, 1,12,4) = [  0,  0,  0]
 tetra_shifts(:, 2,12,4) = [ -1, -1,  1]
 tetra_shifts(:, 3,12,4) = [ -1,  0,  1]
 tetra_shifts(:, 4,12,4) = [  0,  0,  1]
 tetra_shifts(:, 1,13,4) = [  0,  0,  0]
 tetra_shifts(:, 2,13,4) = [  0,  0, -1]
 tetra_shifts(:, 3,13,4) = [  1,  0, -1]
 tetra_shifts(:, 4,13,4) = [  1,  1, -1]
 tetra_shifts(:, 1,14,4) = [  0,  0,  0]
 tetra_shifts(:, 2,14,4) = [  0,  0, -1]
 tetra_shifts(:, 3,14,4) = [  0,  1, -1]
 tetra_shifts(:, 4,14,4) = [  1,  1, -1]
 tetra_shifts(:, 1,15,4) = [  0,  0,  0]
 tetra_shifts(:, 2,15,4) = [  1,  0,  0]
 tetra_shifts(:, 3,15,4) = [  1,  0, -1]
 tetra_shifts(:, 4,15,4) = [  1,  1, -1]
 tetra_shifts(:, 1,16,4) = [  0,  0,  0]
 tetra_shifts(:, 2,16,4) = [  0,  1,  0]
 tetra_shifts(:, 3,16,4) = [  0,  1, -1]
 tetra_shifts(:, 4,16,4) = [  1,  1, -1]
 tetra_shifts(:, 1,17,4) = [  0,  0,  0]
 tetra_shifts(:, 2,17,4) = [  1,  0,  0]
 tetra_shifts(:, 3,17,4) = [  1,  1,  0]
 tetra_shifts(:, 4,17,4) = [  1,  1, -1]
 tetra_shifts(:, 1,18,4) = [  0,  0,  0]
 tetra_shifts(:, 2,18,4) = [  0,  1,  0]
 tetra_shifts(:, 3,18,4) = [  1,  1,  0]
 tetra_shifts(:, 4,18,4) = [  1,  1, -1]
 tetra_shifts(:, 1,19,4) = [  0,  0,  0]
 tetra_shifts(:, 2,19,4) = [  0,  0, -1]
 tetra_shifts(:, 3,19,4) = [  0,  1, -1]
 tetra_shifts(:, 4,19,4) = [ -1,  0,  0]
 tetra_shifts(:, 1,20,4) = [  0,  0,  0]
 tetra_shifts(:, 2,20,4) = [  0,  1,  0]
 tetra_shifts(:, 3,20,4) = [  0,  1, -1]
 tetra_shifts(:, 4,20,4) = [ -1,  0,  0]
 tetra_shifts(:, 1,21,4) = [  0,  0,  0]
 tetra_shifts(:, 2,21,4) = [  0,  0, -1]
 tetra_shifts(:, 3,21,4) = [  1,  0, -1]
 tetra_shifts(:, 4,21,4) = [  0, -1,  0]
 tetra_shifts(:, 1,22,4) = [  0,  0,  0]
 tetra_shifts(:, 2,22,4) = [  1,  0,  0]
 tetra_shifts(:, 3,22,4) = [  1,  0, -1]
 tetra_shifts(:, 4,22,4) = [  0, -1,  0]
 tetra_shifts(:, 1,23,4) = [  0,  0,  0]
 tetra_shifts(:, 2,23,4) = [  0,  0, -1]
 tetra_shifts(:, 3,23,4) = [ -1, -1,  0]
 tetra_shifts(:, 4,23,4) = [  0, -1,  0]
 tetra_shifts(:, 1,24,4) = [  0,  0,  0]
 tetra_shifts(:, 2,24,4) = [  0,  0, -1]
 tetra_shifts(:, 3,24,4) = [ -1, -1,  0]
 tetra_shifts(:, 4,24,4) = [ -1,  0,  0]

 ! These shifts are taken from previous tetrahedron implmentation MJV and BXU
 ! TODO: implement shifts for the other diagonals
 tetra_shifts_6(:,1,1,1) = [0,0,0]
 tetra_shifts_6(:,2,1,1) = [1,0,0]
 tetra_shifts_6(:,3,1,1) = [0,1,0]
 tetra_shifts_6(:,4,1,1) = [1,0,1]
 tetra_shifts_6(:,1,2,1) = [1,0,0]
 tetra_shifts_6(:,2,2,1) = [1,1,0]
 tetra_shifts_6(:,3,2,1) = [0,1,0]
 tetra_shifts_6(:,4,2,1) = [1,0,1]
 tetra_shifts_6(:,1,3,1) = [0,1,0]
 tetra_shifts_6(:,2,3,1) = [1,1,0]
 tetra_shifts_6(:,3,3,1) = [1,0,1]
 tetra_shifts_6(:,4,3,1) = [1,1,1]
 tetra_shifts_6(:,1,4,1) = [0,0,0]
 tetra_shifts_6(:,2,4,1) = [0,1,0]
 tetra_shifts_6(:,3,4,1) = [0,0,1]
 tetra_shifts_6(:,4,4,1) = [1,0,1]
 tetra_shifts_6(:,1,5,1) = [0,0,1]
 tetra_shifts_6(:,2,5,1) = [1,0,1]
 tetra_shifts_6(:,3,5,1) = [0,1,0]
 tetra_shifts_6(:,4,5,1) = [0,1,1]
 tetra_shifts_6(:,1,6,1) = [0,1,0]
 tetra_shifts_6(:,2,6,1) = [1,0,1]
 tetra_shifts_6(:,3,6,1) = [0,1,1]
 tetra_shifts_6(:,4,6,1) = [1,1,1]

 main_diagonals(:,1) = [ 1, 1, 1] ! 0-7
 main_diagonals(:,2) = [-1, 1, 1] ! 1-6
 main_diagonals(:,3) = [ 1,-1, 1] ! 2-5
 main_diagonals(:,4) = [ 1, 1,-1] ! 3-4

 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 tetra%nkibz = nkpt_ibz
 tetra%nkbz = nkpt_fullbz
 ierr = 0
 tetra_count = 0

 ! Determine the smallest diagonal in k-space
 min_length = huge(min_length)
 do idiag = 1,4
   diag(:) = gprimd(:,1)*main_diagonals(1,idiag)+&
             gprimd(:,2)*main_diagonals(2,idiag)+&
             gprimd(:,3)*main_diagonals(3,idiag)
   length = sqrt(diag(1)*diag(1) + diag(2)*diag(2) + diag(3)*diag(3))
   if (length < min_length) then
     min_length = length
     min_idiag = idiag
   end if
 end do

 ! HM TODO: Avoid mkkptrank and map the k-point grid to indexes
 ! Make full k-point rank arrays
 call mkkptrank(kpt_fullbz,nkpt_fullbz,kptrank_t)

 !
 ! HM (13/04/2019): I implement two different versions:
 ! 1. I only use 24 tetrahedra around the IBZ k-point
 ! following the approach of A. Togo (phonopy, spglib, kspclib).
 ! 2. I generate tetrahedra on the full Brillouin zone
 ! and keep track of how many are contributing to the IBZ point and the multiplicities,
 ! these can be more than 24 tetrahedra.
 ! This is equivalent to what Matthieu implemented but avoids large memory allocations.
 !
 ! The two implementations differ specially when using low k-point sampling
 ! (the second yields the same results using IBZ or FBZ, the first one not).
 ! For large sampling the two approaches yield similar results, with the first
 ! one using less memory, faster to generate and compute
 !
#if 0
 ! For each k-point in the IBZ store 24 tetrahedra each refering to 4 k-points
 ABI_MALLOC(tetra%unique_tetra,(tetra%nkibz))
 do ihash=1,tetra%nkibz
   ABI_MALLOC(tetra%unique_tetra(ihash)%indexes,(0:4,24))
   tetra%unique_tetra(ihash)%indexes = 0
 end do
 ! For each k-point in the IBZ
 do ikibz=1,tetra%nkibz
   !if (mod(ikibz,nprocs) /= my_rank) cycle
   k1 = kpt_ibz(:,ikibz)
   tetra_loop: do itetra=1,24
     do isummit=1,4
       ! Find the index of the neighbouring k-points in the BZ
       k2 = k1 + tetra_shifts(1,isummit,itetra,min_idiag)*klatt(:,1) + &
                 tetra_shifts(2,isummit,itetra,min_idiag)*klatt(:,2) + &
                 tetra_shifts(3,isummit,itetra,min_idiag)*klatt(:,3)
       ! Find full kpoint which is summit isummit of tetrahedron itetra around full kpt ikpt_full !
       call get_rank_1kpt(k2,symrankkpt,kptrank_t)
       ikpt2 = kptrank_t%invrank(symrankkpt)
       ! Find the index of those points in the BZ and IBZ
       tetra_ibz(isummit) = bz2ibz(ikpt2)
     end do
     ! Sort index of irr k-point edges (need this so the comparison works)
     call sort_4tetra_int(tetra_ibz)

     ! Store only unique tetrahedra
     ! Compute a very simple hash for each tetrahedron
     ihash = mod(sum(tetra_ibz),tetra%nkibz)+1
     ! Loop over all tetrahedrons that contain this ikibz as first element
     do jtetra=1,tetra_count(ihash)
       ! if tetrahedron already exists add multiplicity
       if (tetra%unique_tetra(ihash)%indexes(1,jtetra)/=tetra_ibz(1)) cycle
       if (tetra%unique_tetra(ihash)%indexes(2,jtetra)/=tetra_ibz(2)) cycle
       if (tetra%unique_tetra(ihash)%indexes(3,jtetra)/=tetra_ibz(3)) cycle
       if (tetra%unique_tetra(ihash)%indexes(4,jtetra)/=tetra_ibz(4)) cycle
       tetra%unique_tetra(ihash)%indexes(0,jtetra) = tetra%unique_tetra(ihash)%indexes(0,jtetra)+1
       cycle tetra_loop
     end do
     ! Otherwise store new tetrahedron
     tetra_count(ihash) = tetra_count(ihash)+1
     max_ntetra = size(tetra%unique_tetra(ihash)%indexes,2)
     ! The contents don't fit the array so I have to resize it
     if (tetra_count(ihash)>max_ntetra) then
       ABI_MALLOC(indexes,(0:4,max_ntetra+TETRA_STEP))
       indexes(0:4,:max_ntetra) = tetra%unique_tetra(ihash)%indexes
       indexes(:,max_ntetra+1:) = 0
       ABI_FREE(tetra%unique_tetra(ihash)%indexes)
       tetra%unique_tetra(ihash)%indexes => indexes
     end if
     tetra%unique_tetra(ihash)%indexes(1:,tetra_count(ihash)) = tetra_ibz(:)
     tetra%unique_tetra(ihash)%indexes(0, tetra_count(ihash)) = 1
   end do tetra_loop
 end do
#else
 min_idiag = 1
 ! For each k-point in the IBZ store 24 tetrahedra each refering to 4 k-points
 ABI_MALLOC(tetra%unique_tetra,(tetra%nkibz))
 do ikibz=1,tetra%nkibz
   ABI_MALLOC(tetra%unique_tetra(ikibz)%indexes,(0:4,TETRA_SIZE))
   tetra%unique_tetra(ikibz)%indexes=0
 end do
 ! For each k-point in the BZ
 do ikbz=1,tetra%nkbz
   k1 = kpt_fullbz(:,ikbz)
   tetra_loop: do itetra=1,6
     ! Determine tetrahedron
     do isummit=1,4
       ! Find the index of the neighbouring k-points in the BZ
       k2 = k1 + tetra_shifts_6(1,isummit,itetra,min_idiag)*klatt(:,1) + &
                 tetra_shifts_6(2,isummit,itetra,min_idiag)*klatt(:,2) + &
                 tetra_shifts_6(3,isummit,itetra,min_idiag)*klatt(:,3)
       ! Find full kpoint which is summit isummit of tetrahedron itetra around full kpt ikpt_full !
       call get_rank_1kpt(k2,symrankkpt,kptrank_t)
       ikpt2 = kptrank_t%invrank(symrankkpt)
       ! Find the index of those points in the BZ and IBZ
       tetra_ibz(isummit) = bz2ibz(ikpt2)
     end do
     ! Sort index of irr k-point edges (need this so the comparison works)
     call sort_4tetra_int(tetra_ibz)

     ! Store only unique tetrahedra
     ! Compute a very simple hash for each tetrahedron
     ihash = mod(sum(tetra_ibz),tetra%nkibz)+1
     ! Loop over all tetrahedrons that contain this ikibz as first element
     do jtetra=1,tetra_count(ihash)
       ! if tetrahedron already exists add multiplicity
       if (tetra%unique_tetra(ihash)%indexes(1,jtetra)/=tetra_ibz(1)) cycle
       if (tetra%unique_tetra(ihash)%indexes(2,jtetra)/=tetra_ibz(2)) cycle
       if (tetra%unique_tetra(ihash)%indexes(3,jtetra)/=tetra_ibz(3)) cycle
       if (tetra%unique_tetra(ihash)%indexes(4,jtetra)/=tetra_ibz(4)) cycle
       tetra%unique_tetra(ihash)%indexes(0,jtetra) = tetra%unique_tetra(ihash)%indexes(0,jtetra)+1
       cycle tetra_loop
     end do
     ! Otherwise store new tetrahedron
     tetra_count(ihash) = tetra_count(ihash)+1
     max_ntetra = size(tetra%unique_tetra(ihash)%indexes,2)
     ! The contents don't fit the array so I have to resize it
     if (tetra_count(ihash)>max_ntetra) then
       ABI_MALLOC(indexes,(0:4,max_ntetra+TETRA_STEP))
       indexes(0:4,:max_ntetra) = tetra%unique_tetra(ihash)%indexes
       indexes(:,max_ntetra+1:) = 0
       ABI_FREE(tetra%unique_tetra(ihash)%indexes)
       tetra%unique_tetra(ihash)%indexes => indexes
     end if
     tetra%unique_tetra(ihash)%indexes(1:,tetra_count(ihash)) = tetra_ibz(:)
     tetra%unique_tetra(ihash)%indexes(0, tetra_count(ihash)) = 1
   end do tetra_loop
 end do
#endif
 call destroy_kptrank(kptrank_t)

 ! Do some maintenance: free unused memory and count tetrahedra per IBZ point
 tetra_count = 0
 do ihash=1,tetra%nkibz
   ntetra = count(tetra%unique_tetra(ihash)%indexes(0,:)>0)
   ! Allocate array with right size
   ABI_MALLOC(indexes,(0:4,ntetra))
   indexes = tetra%unique_tetra(ihash)%indexes(:,:ntetra)
   ABI_FREE(tetra%unique_tetra(ihash)%indexes)
   tetra%unique_tetra(ihash)%indexes => indexes
   ! Count number of tetrahedra per IBZ point
   do itetra=1,ntetra
     tetra_mibz = tetra%unique_tetra(ihash)%indexes(:,itetra)
     do isummit=1,4
       ikibz = tetra_mibz(isummit)
       tetra_count(ikibz) = tetra_count(ikibz) + 1
     end do
   end do
 end do

 ! Allocate IBZ to tetrahedron mapping
 ABI_MALLOC(tetra%ibz,(tetra%nkibz))
 do ikibz=1,tetra%nkibz
   tetra%ibz(ikibz)%tetra_count = tetra_count(ikibz)
   ABI_MALLOC(tetra%ibz(ikibz)%tetra,(tetra_count(ikibz)))
 end do

 ! Create mapping from IBZ to unique tetrahedra
 tetra_count = 0
 do ihash=1,tetra%nkibz
   ntetra = size(tetra%unique_tetra(ihash)%indexes,2)
   do itetra=1,ntetra
     tetra_mibz = tetra%unique_tetra(ihash)%indexes(:,itetra)
     do isummit=1,4
       ikibz = tetra_mibz(isummit)
       tetra_count(ikibz) = tetra_count(ikibz) + 1
       tetra%ibz(ikibz)%tetra(tetra_count(ikibz))%p(0:4) => tetra%unique_tetra(ihash)%indexes(:,itetra)
     end do
   end do
 end do

 ! Sum the multiplicity
 do ikibz=1,tetra%nkibz
   tetra_total = 0
   do itetra=1,tetra%ibz(ikibz)%tetra_count
     tetra_total = tetra_total + tetra%ibz(ikibz)%tetra(itetra)%p(0)
   end do
   tetra%ibz(ikibz)%tetra_total = tetra_total
 end do

 ! Count unique tetra
 total_ntetra = 0
 do ihash=1,tetra%nkibz
   ntetra = size(tetra%unique_tetra(ihash)%indexes,2)
   total_ntetra = total_ntetra + ntetra
 end do
 tetra%nunique_tetra = total_ntetra

 ! Count IBZ tetra
 total_ntetra = 0
 do ikibz=1,tetra%nkibz
   ntetra = tetra%ibz(ikibz)%tetra_count
   total_ntetra = total_ntetra + ntetra
 end do
 tetra%nibz_tetra = total_ntetra

 ! Compute the weights
 ABI_CALLOC(tetra%ibz_weights,(tetra%nkibz))
 do ikbz=1,nkpt_fullbz
   ikibz = bz2ibz(ikbz)
   tetra%ibz_weights(ikibz) = tetra%ibz_weights(ikibz) + 1
 end do
 tetra%ibz_weights = tetra%ibz_weights / nkpt_fullbz

 ! HM TODO: Need to check where this will be used!
 ! Calculate the volume of the tetrahedra
 rcvol = abs(gprimd(1,1)*(gprimd(2,2)*gprimd(3,3)-gprimd(3,2)*gprimd(2,3))- &
             gprimd(2,1)*(gprimd(1,2)*gprimd(3,3)-gprimd(3,2)*gprimd(1,3))+ &
             gprimd(3,1)*(gprimd(1,2)*gprimd(2,3)-gprimd(2,2)*gprimd(1,3)))

 ! Volume of all tetrahedra should be the same as that of tetra 1
 ! this is the volume of 1 tetrahedron, should be coherent with notation in Lehmann & Taut
 k1(:) = gprimd(:,1)*klatt(1,1) +  gprimd(:,2)*klatt(2,1) +  gprimd(:,3)*klatt(3,1)
 k2(:) = gprimd(:,1)*klatt(1,2) +  gprimd(:,2)*klatt(2,2) +  gprimd(:,3)*klatt(3,2)
 k3(:) = gprimd(:,1)*klatt(1,3) +  gprimd(:,2)*klatt(2,3) +  gprimd(:,3)*klatt(3,3)
 tetra%vv = abs(k1(1)*(k2(2)*k3(3)-k2(3)*k3(2))- &
                k1(2)*(k2(1)*k3(3)-k2(3)*k3(1))+ &
                k1(3)*(k2(1)*k3(2)-k2(2)*k3(1))) / 6.d0 / rcvol

end subroutine htetra_init
!!***

!----------------------------------------------------------------------

!!****f* m_htetrahedron/htetra_print
!! NAME
!! htetra_print
!!
!! FUNCTION
!! write information about the tetrahedra object
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine htetra_print(self)

 type(t_htetrahedron), intent(in) :: self
 integer  :: pointer_size_bits
 real(dp) :: total_size, unique_tetra_size, ibz_pointer_size

 unique_tetra_size = self%nunique_tetra*5*four/1024/1024
 pointer_size_bits = storage_size(self%ibz(1)%tetra(1)%p)
 ibz_pointer_size  = one*self%nibz_tetra*pointer_size_bits/8/1024/1024
 total_size        = ibz_pointer_size+unique_tetra_size
 write(std_out,'(a,i12)')     'unique_tetra      ', self%nunique_tetra
 write(std_out,'(a,f12.1,a)') 'unique_tetra_size ', unique_tetra_size, ' [Mb]'
 write(std_out,'(a,i12)')     'ibz_tetra         ', self%nibz_tetra
 write(std_out,'(a,i12,a)')   'pointer_size      ', pointer_size_bits, ' bits'
 write(std_out,'(a,f12.1,a)') 'ibz_pointer_size  ', ibz_pointer_size, ' [Mb]'
 write(std_out,'(a,f12.1,a)') 'total size        ', total_size, ' [Mb]'

end subroutine htetra_print
!!***

!----------------------------------------------------------------------

!!****f* m_htetrahedron/htetra_free
!! NAME
!! htetra_free
!!
!! FUNCTION
!! deallocate tetrahedra pointers if needed
!!
!! PARENTS
!!      ep_el_weights,ep_fs_weights,ep_ph_weights,gstate,m_ebands,m_epjdos
!!      m_fstab,m_gruneisen,m_phgamma,m_phonons,thmeig,wfk_analyze
!!
!! CHILDREN
!!
!! SOURCE

subroutine htetra_free(tetra)

 type(t_htetrahedron), intent(inout) :: tetra
 integer :: ikibz

 do ikibz=1,tetra%nkibz
   ABI_SFREE(tetra%ibz(ikibz)%tetra)
   ABI_FREE(tetra%unique_tetra(ikibz)%indexes)
 end do
 ABI_SFREE(tetra%ibz)
 ABI_SFREE(tetra%unique_tetra)
 ABI_SFREE(tetra%ibz_weights)

end subroutine htetra_free
!!***

!----------------------------------------------------------------------

!!****f* m_htetrahedron/get_onetetra_
!! NAME
!! get_onetetra_
!!
!! FUNCTION
!! Private function to calculate the contributions to the weights due to a single tetrahedron.
!! Extracted from get_tetra_weight
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

pure subroutine get_onetetra_(tetra,eigen_1tetra,energies,nene,max_occ,bcorr, &
                              tweight,dweight)

!Arguments ------------------------------------
!scalars
 integer,intent(in)   :: nene,bcorr
 real(dp) ,intent(in) :: max_occ
 real(dp) ,intent(in) :: energies(nene)
 type(t_htetrahedron), intent(in) :: tetra
!arrays
 real(dp), intent(out) :: tweight(4, nene)
 real(dp), intent(out) :: dweight(4, nene)
 real(dp),intent(in)  :: eigen_1tetra(4)

!Local variables-------------------------------
 integer :: ieps
 real(dp) :: cc,cc1,cc2,cc3
 real(dp) :: dcc1de,dcc2de,dcc3de,dccde,eps
 real(dp) :: e21,e31,e32,e41,e42,e43
 real(dp) :: inv_e32,inv_e41,inv_e42,inv_e43,inv_e21,inv_e31
 real(dp) :: deleps1,deleps2,deleps3,deleps4
 real(dp) :: e1,e2,e3,e4
 real(dp) :: invepsum, cc_pre, dccde_pre
 real(dp) :: cc1_pre, cc2_pre, cc3_pre
 real(dp) :: dccde_tmp
 real(dp) :: bcorr_fact,volconst

! *********************************************************************

 ! Factor of 1/6 from the volume of the tetrahedron
 ! The rest is accounted for from the kpoint weights
 volconst = max_occ!/4.d0/6.d0

 ! This is output
 tweight = zero; dweight = zero

 ! all notations are from Blochl PRB 49 16223 [[cite:Bloechl1994a]] Appendix B
 e1 = eigen_1tetra(1)
 e2 = eigen_1tetra(2)
 e3 = eigen_1tetra(3)
 e4 = eigen_1tetra(4)
 e21 = e2-e1
 e31 = e3-e1
 e41 = e4-e1
 e32 = e3-e2
 e42 = e4-e2
 e43 = e4-e3
 inv_e21 = zero; if (e21 > tol14) inv_e21 = 1.d0 / e21
 inv_e31 = zero; if (e31 > tol14) inv_e31 = 1.d0 / e31
 inv_e41 = zero; if (e41 > tol14) inv_e41 = 1.d0 / e41
 inv_e32 = zero; if (e32 > tol14) inv_e32 = 1.d0 / e32
 inv_e42 = zero; if (e42 > tol14) inv_e42 = 1.d0 / e42
 inv_e43 = zero; if (e43 > tol14) inv_e43 = 1.d0 / e43

 do ieps=1,nene
   eps = energies(ieps)

   !
   ! eps < e1 nothing to do
   !
   if (eps < e1) cycle

   !
   ! e1 < eps < e2
   !
   if (eps < e2) then
     deleps1 = eps-e1
     invepsum = inv_e21+inv_e31+inv_e41

     ! Heaviside
     cc = volconst*inv_e21*inv_e31*inv_e41*deleps1**3
     tweight(1,ieps) = cc*(4.d0-deleps1*invepsum)
     tweight(2,ieps) = cc*deleps1*inv_e21
     tweight(3,ieps) = cc*deleps1*inv_e31
     tweight(4,ieps) = cc*deleps1*inv_e41

     ! Delta
     dccde_pre = 3.d0*volconst*inv_e21*inv_e31*inv_e41
     dccde = dccde_pre*deleps1**2
     dweight(1,ieps) = dccde*(4.d0-deleps1*invepsum)-cc*invepsum
     dweight(2,ieps) = (dccde*deleps1+cc) * inv_e21
     dweight(3,ieps) = (dccde*deleps1+cc) * inv_e31
     dweight(4,ieps) = (dccde*deleps1+cc) * inv_e41

     if (bcorr == 1) then
       ! bxu, correction terms based on Bloechl's paper
       bcorr_fact = 4.d0/40.d0*dccde_pre*deleps1*deleps1
       tweight(1,ieps) = tweight(1,ieps) + bcorr_fact*( e21+e31+e41)
       tweight(2,ieps) = tweight(2,ieps) + bcorr_fact*(-e21+e32+e42)
       tweight(3,ieps) = tweight(3,ieps) + bcorr_fact*(-e31-e32+e43)
       tweight(4,ieps) = tweight(4,ieps) + bcorr_fact*(-e41-e42-e43)

       bcorr_fact = 8.d0/40.d0*dccde_pre*deleps1
       dweight(1,ieps) = dweight(1,ieps) + bcorr_fact*( e21+e31+e41)
       dweight(2,ieps) = dweight(2,ieps) + bcorr_fact*(-e21+e32+e42)
       dweight(3,ieps) = dweight(3,ieps) + bcorr_fact*(-e31-e32+e43)
       dweight(4,ieps) = dweight(4,ieps) + bcorr_fact*(-e41-e42-e43)
     end if

     cycle
   endif

   !
   ! e2 < eps < e3
   !
   if (eps < e3) then
     deleps1 = eps-e1
     deleps2 = eps-e2
     deleps3 = e3-eps
     deleps4 = e4-eps
     cc1_pre = volconst*inv_e31*inv_e41
     cc2_pre = volconst*inv_e41*inv_e32*inv_e31
     cc3_pre = volconst*inv_e42*inv_e32*inv_e41

     ! Heaviside
     cc1 = cc1_pre*deleps1*deleps1
     cc2 = cc2_pre*deleps1*deleps2*deleps3
     cc3 = cc3_pre*deleps2*deleps2*deleps4

     tweight(1,ieps) = (cc1)+&
                       (cc1+cc2)*deleps3*inv_e31+&
                       (cc1+cc2+cc3)*deleps4*inv_e41
     tweight(2,ieps) = (cc1+cc2+cc3)+&
                       (cc2+cc3)*deleps3*inv_e32+&
                           (cc3)*deleps4*inv_e42
     tweight(3,ieps) = (cc1+cc2)*deleps1*inv_e31+&
                       (cc2+cc3)*deleps2*inv_e32
     tweight(4,ieps) = (cc1+cc2+cc3)*deleps1*inv_e41+&
                                   (cc3)*deleps2*inv_e42

     ! Delta
     dcc1de = 2.d0*cc1_pre*(     deleps1)
     dcc2de =      cc2_pre*(    -deleps1*deleps2+deleps1*deleps3+deleps2*deleps3)
     dcc3de =      cc3_pre*(2.d0*deleps2*deleps4-deleps2*deleps2)
     dweight(1,ieps) = dcc1de+&
                       ((dcc1de+dcc2de)*deleps3-(cc1+cc2))*inv_e31+&
                       ((dcc1de+dcc2de+dcc3de)*deleps4-(cc1+cc2+cc3))*inv_e41
     dweight(2,ieps) = (dcc1de+dcc2de+dcc3de)+&
                       ((dcc2de+dcc3de)*deleps3-(cc2+cc3))*inv_e32+&
                                      (dcc3de*deleps4-cc3)*inv_e42
     dweight(3,ieps) = ((dcc1de+dcc2de)*deleps1+(cc1+cc2))*inv_e31+&
                       ((dcc2de+dcc3de)*deleps2+(cc2+cc3))*inv_e32
     dweight(4,ieps) = ((dcc1de+dcc2de+dcc3de)*deleps1+(cc1+cc2+cc3))*inv_e41+&
                                                        (dcc3de*deleps2+cc3)*inv_e42

     if (bcorr == 1) then
       ! bxu, correction terms based on Bloechl's paper
       ! The correction terms may cause the dweight become negative
       bcorr_fact = 4.d0/40.d0*cc1_pre*(3.d0*e21+6.d0*deleps2-3.d0*(e31+e42)*deleps2*deleps2*inv_e32*inv_e42)
       tweight(1,ieps) = tweight(1,ieps) + bcorr_fact*( e21+e31+e41)
       tweight(2,ieps) = tweight(2,ieps) + bcorr_fact*(-e21+e32+e42)
       tweight(3,ieps) = tweight(3,ieps) + bcorr_fact*(-e31-e32+e43)
       tweight(4,ieps) = tweight(4,ieps) + bcorr_fact*(-e41-e42-e43)

       bcorr_fact = 4.d0/40.d0*cc1_pre*(6.d0-6.d0*(e31+e42)*deleps2*inv_e32*inv_e42)
       dweight(1,ieps) = dweight(1,ieps) + bcorr_fact*( e21+e31+e41)
       dweight(2,ieps) = dweight(2,ieps) + bcorr_fact*(-e21+e32+e42)
       dweight(3,ieps) = dweight(3,ieps) + bcorr_fact*(-e31-e32+e43)
       dweight(4,ieps) = dweight(4,ieps) + bcorr_fact*(-e41-e42-e43)
     end if

     cycle
   endif


   !
   ! e3 < eps < e4
   !
   if (eps < e4) then
     deleps4 = e4-eps
     invepsum = inv_e41+inv_e42+inv_e43

     ! Heaviside
     cc_pre = volconst*inv_e41*inv_e42*inv_e43
     cc = cc_pre*deleps4**3
     tweight(1,ieps) = volconst - deleps4*cc*inv_e41
     tweight(2,ieps) = volconst - deleps4*cc*inv_e42
     tweight(3,ieps) = volconst - deleps4*cc*inv_e43
     tweight(4,ieps) = volconst - cc*(4.d0-deleps4*invepsum)

     ! Delta
     dccde = -3.d0*cc_pre*deleps4**2
     dccde_tmp = dccde*deleps4 + cc
     dweight(1,ieps) = -dccde_tmp * inv_e41
     dweight(2,ieps) = -dccde_tmp * inv_e42
     dweight(3,ieps) = -dccde_tmp * inv_e43
     dweight(4,ieps) = -4.d0*dccde + dccde_tmp*invepsum

     if (bcorr == 1) then
       ! bxu, correction terms based on Bloechl's paper
       ! The correction terms may cause the dweight become negative
       bcorr_fact = 12.d0/40.d0*cc_pre*deleps4*deleps4
       tweight(1,ieps) = tweight(1,ieps) + bcorr_fact*( e21+e31+e41)
       tweight(2,ieps) = tweight(2,ieps) + bcorr_fact*(-e21+e32+e42)
       tweight(3,ieps) = tweight(3,ieps) + bcorr_fact*(-e31-e32+e43)
       tweight(4,ieps) = tweight(4,ieps) + bcorr_fact*(-e41-e42-e43)

       bcorr_fact = - 24.d0/40.d0*cc_pre*deleps4
       dweight(1,ieps) = dweight(1,ieps) + bcorr_fact*( e21+e31+e41)
       dweight(2,ieps) = dweight(2,ieps) + bcorr_fact*(-e21+e32+e42)
       dweight(3,ieps) = dweight(3,ieps) + bcorr_fact*(-e31-e32+e43)
       dweight(4,ieps) = dweight(4,ieps) + bcorr_fact*(-e41-e42-e43)
     end if

     cycle
   endif

   !
   ! e4 < eps
   !
   if (e4 < eps) then

     ! Heaviside
     tweight(:,ieps:) = volconst

     ! Delta unchanged by this tetrahedron
     exit
   end if

   !
   !  if we have a fully degenerate tetrahedron,
   !  1) the tweight is a Heaviside (step) function, which is correct above, but
   !  2) the dweight should contain a Dirac function
   !
 end do

end subroutine get_onetetra_
!!***

!----------------------------------------------------------------------

!!****f* m_htetrahedron/htetra_get_onewk_wvals
!! NAME
!! htetra_get_onewk_wvals
!!
!! FUNCTION
!! Calculate integration weights and their derivatives for a single k-point in the IBZ.
!!
!! INPUTS
!! tetra<t_htetrahedron>=Object with tables for tetrahedron method.
!! ik_ibz=Index of the k-point in the IBZ array
!! bcorr=1 to include Blochl correction else 0.
!! nw=number of energies in wvals
!! nibz=number of irreducible kpoints
!! wvals(nw)=Frequency points.
!! eigen_ibz(nkibz)=eigenenergies for each k point
!!
!! OUTPUT
!!  weights(nw,2) = integration weights for
!!    Dirac delta (derivative of theta wrt energy) and Theta (Heaviside function)
!!    for a given (band, k-point, spin).
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine htetra_get_onewk_wvals(tetra, ik_ibz, bcorr, nw, wvals, max_occ, nkibz, eig_ibz, weights)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_ibz,nw,nkibz,bcorr
 real(dp) ,intent(in) :: max_occ
 type(t_htetrahedron), intent(in) :: tetra
!arrays
 real(dp),intent(in) :: wvals(nw)
 real(dp),intent(in) :: eig_ibz(nkibz)
 real(dp),intent(out) :: weights(nw, 2)

!Local variables-------------------------------
!scalars
 real(dp) :: tweight
 integer  :: itetra,isummit,tetra_count,tetra_total
!arrays
 integer  :: ind_ibz(4)
 real(dp) :: eigen_1tetra(4)
 real(dp) :: tweight_tmp(4,nw),dtweightde_tmp(4,nw)

! *********************************************************************

 weights = zero

 ! For each tetrahedron that belongs to this k-point
 tetra_count = tetra%ibz(ik_ibz)%tetra_count
 tetra_total = tetra%ibz(ik_ibz)%tetra_total
 do itetra=1,tetra_count

   tweight = one*tetra%ibz(ik_ibz)%tetra(itetra)%p(0)/tetra_total
   do isummit=1,4
     ! Get mapping of each summit to eig_ibz
     ind_ibz(isummit) = tetra%ibz(ik_ibz)%tetra(itetra)%p(isummit)
     eigen_1tetra(isummit) = eig_ibz(ind_ibz(isummit))
   end do

   ! Sort energies before calling get_onetetra_
   !call sort_tetra(4, eigen_1tetra, ind_ibz, tol14)
   call sort_4tetra(eigen_1tetra, ind_ibz)
   ! HM: Here we should only compute what we will use!
   call get_onetetra_(tetra, eigen_1tetra, wvals, nw, max_occ, bcorr, tweight_tmp, dtweightde_tmp)

   ! Accumulate contributions to ik_ibz (there might be multiple vertexes that map onto ik_ibz)
   do isummit=1,4
     if (ind_ibz(isummit) /= ik_ibz) cycle
     weights(:,1) = weights(:,1) + dtweightde_tmp(isummit,:)*tweight
     weights(:,2) = weights(:,2) + tweight_tmp(isummit,:)   *tweight
     ! HM: This exit is important, avoids summing the same contribution more than once
     exit
   end do
 end do ! itetra

end subroutine htetra_get_onewk_wvals
!!***

!----------------------------------------------------------------------

!!****f* m_htetrahedron/tetra_get_onewk
!! NAME
!! tetra_get_onewk
!!
!! FUNCTION
!! Calculate integration weights and their derivatives for a single k-point in the IBZ.
!! Same as above but different calling arguments.
!! HM: The above is prefered but I keep this one to ease the transition
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

subroutine htetra_get_onewk(tetra,ik_ibz,bcorr,nw,nkibz,eig_ibz,&
                           enemin,enemax,max_occ,weights)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_ibz,nw,nkibz,bcorr
 type(t_htetrahedron), intent(in) :: tetra
 real(dp) ,intent(in) :: enemin,enemax,max_occ
!arrays
 real(dp),intent(in) :: eig_ibz(nkibz)
 real(dp),intent(out) :: weights(nw,2)

!Local variables-------------------------------
!scalars
 integer :: itetra,isummit
!arrays
 integer :: ind_ibz(4)
 real(dp) :: wvals(nw)
 real(dp) :: eigen_1tetra(4)
 real(dp) :: tweight_tmp(nw,4),dtweightde_tmp(nw,4)

! *********************************************************************

 weights = zero
 wvals = linspace(enemin,enemax,nw)
 call htetra_get_onewk_wvals(tetra, ik_ibz, bcorr, nw, wvals, max_occ, nkibz, eig_ibz, weights)

end subroutine htetra_get_onewk
!!***

!----------------------------------------------------------------------

!!****f* m_htetrahedron/htetra_get_onewk_wvals_zinv
!! NAME
!! htetra_get_onewk_wvals_zinv
!!
!! FUNCTION
!! Calculate integration weights for 1/(z-E(k)) for a single k-point in the IBZ.
!! Using either the implementation from:
!! S. Kaprzyk, Computer Physics Communications 183, 347 (2012).
!! or (TODO)
!! P. Lambin and J.P. Vigneron, Phys. Rev. B 29, 3430 (1984).
!!
!! INPUTS
!! tetra<t_htetrahedron>=Object with tables for tetrahedron method.
!! ik_ibz=Index of the k-point in the IBZ array
!! bcorr=1 to include Blochl correction else 0.
!! nw=number of energies in wvals
!! nibz=number of irreducible kpoints
!! wvals(nw)=Frequency points.
!! eigen_ibz(nkibz)=eigenenergies for each k point
!!
!! OUTPUT
!!  weights(nw,2) = integration weights for
!!    Dirac delta (derivative of theta wrt energy) and Theta (Heaviside function)
!!    for a given (band, k-point, spin).
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine htetra_get_onewk_wvals_zinv(tetra, ik_ibz, nz, zvals, max_occ, nkibz, eig_ibz, cweights)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_ibz,nz,nkibz
 real(dp) ,intent(in) :: max_occ
 type(t_htetrahedron), intent(in) :: tetra
!arrays
 complex(dp),intent(in) :: zvals(nz)
 real(dp),intent(in) :: eig_ibz(nkibz)
 complex(dp),intent(out) :: cweights(nz, 2)

!Local variables-------------------------------
!scalars
 integer  :: itetra,isummit,tetra_total,tetra_count,iz
 real(dp) :: tweight
!arrays
 integer  :: ind_ibz(4)
 real(dp) :: eigen_1tetra(4)
 complex(dp) :: verm(4), verl(4), verli(4)
! *********************************************************************

 cweights = zero

 ! TODO: check multiplicity here!
 ! For each tetrahedron that belongs to this k-point
 tetra_count = tetra%ibz(ik_ibz)%tetra_count
 tetra_total = tetra%ibz(ik_ibz)%tetra_total
 do itetra=1,24

   tweight = one*tetra%ibz(ik_ibz)%tetra(itetra)%p(0)/tetra_total
   do isummit=1,4
     ! Get mapping of each summit to eig_ibz
     ind_ibz(isummit) = tetra%ibz(ik_ibz)%tetra(itetra)%p(isummit)
     eigen_1tetra(isummit) = eig_ibz(ind_ibz(isummit))
   end do

   ! Loop over frequencies
   do iz=1,nz
     verm = zvals(iz) - eigen_1tetra
     call SIM0TWOI(VERL, VERLI, VERM)
     do isummit=1,4
       if (ind_ibz(isummit) /= ik_ibz) cycle
       cweights(iz,1) = cweights(iz,1) + verl(isummit)* tweight
       cweights(iz,2) = cweights(iz,2) + verli(isummit)*tweight
       ! HM: This exit is important, avoids summing the same contribution more than once
       exit
     end do
   end do
 end do ! itetra

end subroutine htetra_get_onewk_wvals_zinv
!!***

!----------------------------------------------------------------------

!!****f* m_htetrahedron/htetra_blochl_weights
!! NAME
!!  htetra_blochl_weights
!!
!! FUNCTION
!!   Emulates the behaviour of the previous tetrahedron implementation.
!!   HM: I find that in many routines its better to change the implementation
!!   and accumulate the tetrahedron weights in the same way as the
!!   gaussian smearing weights using htetra_get_onewk_wvals. However this requires
!!   some refactoring of the code. I provide this routine to make it easier
!!   to transition to the new tetrahedron implementation without refactoring.
!!   New implementations should try to use htetra_get_onewk_wvals.
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

subroutine htetra_blochl_weights(tetra,eig_ibz,enemin,enemax,max_occ,nw,nkpt,&
  bcorr,tweight,dweight,comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nw,nkpt,bcorr,comm
 type(t_htetrahedron), intent(in) :: tetra
 real(dp) ,intent(in) :: enemax,enemin,max_occ
!arrays
 real(dp) ,intent(in) :: eig_ibz(nkpt)
 real(dp) ,intent(out) :: dweight(nw,nkpt),tweight(nw,nkpt)

!Local variables-------------------------------
!scalars
 integer :: ik_ibz,multiplicity,nprocs,my_rank,ierr,ii
 integer :: tetra_total, tetra_count, itetra, isummit, ihash
 real :: mweight
!arrays
 integer :: ind_ibz(4)
 real(dp) :: eigen_1tetra(4)
 real(dp) :: wvals(nw)
 real(dp) :: dweight_tmp(4,nw),tweight_tmp(4,nw)

! *********************************************************************

 tweight = zero; dweight = zero
 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 wvals = linspace(enemin,enemax,nw)
 tetra_total = 0

 ! For each bucket of tetrahedra
 do ihash=1,tetra%nkibz
   if (mod(ihash,nprocs) /= my_rank) cycle

   ! For each tetrahedron that belongs to this k-point
   tetra_count = size(tetra%unique_tetra(ihash)%indexes,2)
   do itetra=1,tetra_count

     ! Get mapping of each summit to eig_ibz
     do isummit=1,4
       ind_ibz(isummit) = tetra%unique_tetra(ihash)%indexes(isummit,itetra)
       eigen_1tetra(isummit) = eig_ibz(ind_ibz(isummit))
     end do

     ! Sort energies before calling get_onetetra_
     call sort_4tetra(eigen_1tetra, ind_ibz)

     ! Get tetrahedron weights
     call get_onetetra_(tetra, eigen_1tetra, wvals, nw, max_occ, bcorr, tweight_tmp, dweight_tmp)

     ! Acumulate the contributions
     multiplicity = tetra%unique_tetra(ihash)%indexes(0,itetra)
     tetra_total = tetra_total + multiplicity
     do isummit=1,4
       ik_ibz = ind_ibz(isummit)
       tweight(:,ik_ibz) = tweight(:,ik_ibz) + tweight_tmp(isummit,:)*multiplicity
       dweight(:,ik_ibz) = dweight(:,ik_ibz) + dweight_tmp(isummit,:)*multiplicity
     end do
   end do ! itetra
 end do

 ! Rescale with weights
 do ik_ibz=1,tetra%nkibz
   tweight(:,ik_ibz) = tweight(:,ik_ibz)*tetra%ibz_weights(ik_ibz)/tetra%ibz(ik_ibz)%tetra_total
   dweight(:,ik_ibz) = dweight(:,ik_ibz)*tetra%ibz_weights(ik_ibz)/tetra%ibz(ik_ibz)%tetra_total
 end do

 call xmpi_sum(tweight, comm, ierr)
 call xmpi_sum(dweight, comm, ierr)

end subroutine htetra_blochl_weights
!!***

!----------------------------------------------------------------------

!!****f* m_htetrahedron/sort_4tetra
!! NAME
!!  sort_4tetra
!!
!! FUNCTION
!!  Sort double precision array list(4) into ascending numerical order
!!  while making corresponding rearrangement of the integer
!!  array iperm.
!!
!!  Taken from:
!!  https://stackoverflow.com/questions/6145364/sort-4-number-with-few-comparisons
!!
!! INPUTS
!!  list(4) intent(inout) list of double precision numbers to be sorted
!!  perm(4) intent(inout) iperm(i)=i (very important)
!!
!! OUTPUT
!!  list(4) sorted list
!!  perm(4) index of permutation given the right ascending order
!!
!! PARENTS
!!      m_htetrahedron
!!
!! CHILDREN
!!
!! SOURCE


pure subroutine sort_4tetra(list,perm)

 integer,  intent(inout) :: perm(4)
 real(dp), intent(inout) :: list(4)

!Local variables-------------------------------
 integer :: ia,ib,ic,id
 integer :: ilow1,ilow2,ihigh1,ihigh2
 integer :: ilowest,ihighest
 integer :: imiddle1,imiddle2
 real(dp) :: va,vb,vc,vd
 real(dp) :: vlow1,vlow2,vhigh1,vhigh2
 real(dp) :: vlowest,vhighest
 real(dp) :: vmiddle1,vmiddle2

 va = list(1); ia = perm(1)
 vb = list(2); ib = perm(2)
 vc = list(3); ic = perm(3)
 vd = list(4); id = perm(4)

 if (va < vb) then
     vlow1 = va; vhigh1 = vb
     ilow1 = ia; ihigh1 = ib
 else
     vlow1 = vb; vhigh1 = va
     ilow1 = ib; ihigh1 = ia
 endif

 if (vc < vd) then
     vlow2 = vc; vhigh2 = vd
     ilow2 = ic; ihigh2 = id
 else
     vlow2 = vd; vhigh2 = vc
     ilow2 = id; ihigh2 = ic
 endif

 if (vlow1 < vlow2) then
     vlowest  = vlow1; vmiddle1 = vlow2
     ilowest  = ilow1; imiddle1 = ilow2
 else
     vlowest  = vlow2; vmiddle1 = vlow1
     ilowest  = ilow2; imiddle1 = ilow1
 endif

 if (vhigh1 > vhigh2) then
     vhighest = vhigh1; vmiddle2 = vhigh2
     ihighest = ihigh1; imiddle2 = ihigh2
 else
     vhighest = vhigh2; vmiddle2 = vhigh1
     ihighest = ihigh2; imiddle2 = ihigh1
 endif

 if (vmiddle1 < vmiddle2) then
     list = [vlowest,vmiddle1,vmiddle2,vhighest]
     perm = [ilowest,imiddle1,imiddle2,ihighest]
 else
     list = [vlowest,vmiddle2,vmiddle1,vhighest]
     perm = [ilowest,imiddle2,imiddle1,ihighest]
 endif

end subroutine sort_4tetra
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/sort_4tetra_int
!! NAME
!!  sort_4tetra_int
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

pure subroutine sort_4tetra_int(list)

 integer, intent(inout) :: list(4)

!Local variables-------------------------------
 integer :: va,vb,vc,vd
 integer :: vlow1,vlow2,vhigh1,vhigh2
 integer :: vlowest,vhighest
 integer :: vmiddle1,vmiddle2

 va = list(1)
 vb = list(2)
 vc = list(3)
 vd = list(4)

 if (va < vb) then
     vlow1 = va; vhigh1 = vb
 else
     vlow1 = vb; vhigh1 = va
 endif

 if (vc < vd) then
     vlow2 = vc; vhigh2 = vd
 else
     vlow2 = vd; vhigh2 = vc
 endif

 if (vlow1 < vlow2) then
     vlowest  = vlow1; vmiddle1 = vlow2
 else
     vlowest  = vlow2; vmiddle1 = vlow1
 endif

 if (vhigh1 > vhigh2) then
     vhighest = vhigh1; vmiddle2 = vhigh2
 else
     vhighest = vhigh2; vmiddle2 = vhigh1
 endif

 if (vmiddle1 < vmiddle2) then
     list = [vlowest,vmiddle1,vmiddle2,vhighest]
 else
     list = [vlowest,vmiddle2,vmiddle1,vhighest]
 endif

end subroutine sort_4tetra_int
!!***

!----------------------------------------------------------------------

!!****f* m_numeric_tools/linspace
!! NAME
!!  linspace
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

pure function linspace(start,stop,nn)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nn
 real(dp),intent(in) :: start,stop
 real(dp) :: length
 real(dp) :: linspace(nn)

!Local variables-------------------------------
 integer :: ii
! *********************************************************************

 select case (nn)

 case (1:)
  length = stop-start
  do ii=1,nn
   linspace(ii)=start+length*(ii-1)/(nn-1)
  end do

 case (0)
  RETURN

 end select

end function linspace
!!***

!----------------------------------------------------------------------

end module m_htetrahedron
!!***

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
 use m_numeric_tools,   only : linspace
 use m_simtet,          only : sim0onei, SIM0TWOI
 use m_xmpi

implicit none

private
!!***

integer, parameter :: TETRA_SIZE = 6
integer, parameter :: TETRA_STEP = 6

!!****t* m_htetrahedron/t_htetra_bucket
!! NAME
!! t_htetra_bucket
!!
!! FUNCTION
!! Store a bunch of tetrahedra
!!
!! SOURCE
type :: htetra_bucket

  integer,allocatable :: indexes(:,:)

end type htetra_bucket
!!***

!!****t* m_htetrahedron/t_htetrahedron
!! NAME
!! t_htetrahedron
!!
!! FUNCTION
!! tetrahedron geometry object
!!
!! SOURCE
type, public :: t_htetrahedron

  integer :: opt
  ! Option for the generation of tetrahedra

  integer :: nkibz
  ! Number of points in the irreducible Brillouin zone

  integer :: nkbz
  ! Number of points in the full Brillouin zone

  integer :: nbuckets
  ! Number of buckets for the hash table

  integer :: nunique_tetra
  ! Number of unique tetrahedron

  integer :: nibz_tetra
  ! Number of ibz tetrahedron

  integer,allocatable :: tetra_total(:)
  ! Equivalent tetrahedra per kpoint (number of irred tetra times multiplicity)

  integer,allocatable :: tetra_count(:)
  ! Inequivalent tetrahedra per kpoint (number of irred tetra)

  integer,allocatable :: ibz_multiplicity(:)
  ! Multiplicity of each k-point

  real(dp)  :: vv
  ! volume of the tetrahedra

  real(dp) :: klatt(3, 3)
  ! reciprocal of lattice vectors for full kpoint grid

  type(htetra_bucket),allocatable :: ibz(:)
  ! indexes of the tetrahedra for each k-point

  type(htetra_bucket),allocatable :: unique_tetra(:)
  ! indexes of the unique tetrahedra

  !contains

end type t_htetrahedron
!!***

public :: htetra_init            ! Initialize the object
public :: htetra_free            ! Free memory
public :: htetra_print           ! Print information about tetrahedron object
public :: htetra_get_onewk       ! Calculate integration weights and their derivatives for a single k-point in the IBZ.
public :: htetra_get_onewk_wvals ! Similar to tetra_get_onewk_wvals but receives arbitrary list of frequency points.
public :: htetra_get_onewk_wvals_zinv ! Calculate integration weights for 1/(z-E(k)) for a single k-point in the IBZ.
public :: htetra_weights_wvals_zinv ! Same as above but return the weight on all the kpoints by looping over tetrahedra
public :: htetra_wvals_weights   ! Compute delta and theta on a list of energies for all kpoints
public :: htetra_wvals_weights_delta ! Compute delta on a list of energies for all kpoints
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
!!  options=1.generate 24 tetrahedra per k-point
!!            faster but gives different results depending on the IBZ, small error for large grids
!!          2.generate tetrahedra on the FBZ and map to IBZ
!!            slower but same results for IBZ and FBZ.
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

subroutine htetra_init(tetra, bz2ibz, gprimd, klatt, kpt_fullbz, nkpt_fullbz, kpt_ibz, nkpt_ibz, ierr, errorstring, comm, opt)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt_fullbz, nkpt_ibz, comm
 integer,optional,intent(in) :: opt
 integer,intent(out) :: ierr
 character(len=80),intent(out) :: errorstring
 type(t_htetrahedron),intent(out),target :: tetra
!arrays
 integer,intent(in) :: bz2ibz(nkpt_fullbz)
 real(dp),intent(in) :: gprimd(3,3),klatt(3,3),kpt_fullbz(3,nkpt_fullbz),kpt_ibz(3,nkpt_ibz)

!Local variables-------------------------------
!scalars
 !type(octree_t) :: oct
 type(kptrank_type) :: kptrank
 integer :: ikpt2,isummit,itetra,jtetra
 integer :: ikibz,ikbz,idiag,ihash,min_idiag,my_rank,nprocs
 integer :: max_ntetra, ntetra
 real(dp) :: rcvol,length,min_length
!arrays
 integer,allocatable,target :: indexes(:,:), tetra_hash_count(:)
 integer :: tetra_ibz(4)
 integer :: tetra_shifts(3,4,24,4)  ! 3 dimensions, 4 summits, 24 tetrahedra, 4 main diagonals
 integer :: tetra_shifts_6(3,4,6,1) ! 3 dimensions, 4 summits, 6 tetrahedra, 4 main diagonals
 integer :: main_diagonals(3,4)
 integer :: tetra_mibz(0:4)
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

 ! These shifts are taken from previous tetrahedron implementation by MJV and BXU
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
 ! HM: this value should be important for performance
 ! more buckets means faster queries for unique tetrahedra
 ! but more memory due to the initial size TETRA_SIZE of the buckets
 ! when changing the number of buckets one should also change the hash function
 ! to distribute the tetrahedra in the buckets as uniformly as possible
 ! the simplest hash (not the best!) is:
 ! ihash = mod(sum(tetra_ibz),nbuckts)
 ! the value of sum(tetra_ibz) is between 1 and 4*nkibz so I use nkibz nbuckets
 ! a larger number of buckets should speed up finding the irreducible tetrahedra
 ! but with more memory allocated
 tetra%nbuckets = nkpt_ibz
 ierr = 0
 ABI_CALLOC(tetra_hash_count,(tetra%nbuckets))

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
 !oct = octree_init(kpt_fullbz,2**4,[-one,-one,-one],[two,two,two])
 call mkkptrank(kpt_fullbz,nkpt_fullbz,kptrank)

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
 ABI_MALLOC(tetra%unique_tetra,(tetra%nbuckets))
 do ihash=1,tetra%nbuckets
   ABI_CALLOC(tetra%unique_tetra(ihash)%indexes,(0:4,TETRA_SIZE))
 end do
 tetra%opt = 2; if (present(opt)) tetra%opt = opt
 select case(tetra%opt)
 case(1)
   ! For each k-point in the IBZ store 24 tetrahedra each refering to 4 k-points
   do ikibz=1,tetra%nkibz
     !if (mod(ikibz,nprocs) /= my_rank) cycle
     k1 = kpt_ibz(:,ikibz)
     tetra_loop1: do itetra=1,24
       do isummit=1,4
         ! Find the index of the neighbouring k-points in the BZ
         k2 = k1 + tetra_shifts(1,isummit,itetra,min_idiag)*klatt(:,1) + &
                   tetra_shifts(2,isummit,itetra,min_idiag)*klatt(:,2) + &
                   tetra_shifts(3,isummit,itetra,min_idiag)*klatt(:,3)
         ! Find full kpoint which is summit isummit of tetrahedron itetra around full kpt ikpt_full !
         !ikpt2 = octree_find(oct,k2,dist)
         !ikpt2 = octree_find_nearest_pbc(oct,k2,dist,shift)
         !if (dist>tol12) call exit(1)
         ikpt2 = kptrank_index(kptrank,k2)
         ! Find the index of those points in the BZ and IBZ
         tetra_ibz(isummit) = bz2ibz(ikpt2)
       end do
       ! Sort index of irr k-point edges (need this so the comparison works)
       call sort_4tetra_int(tetra_ibz)

       ! Store only unique tetrahedra
       ! Compute a very simple hash for each tetrahedron
       ihash = compute_hash(tetra,tetra_ibz) !mod(sum(tetra_ibz),tetra%nbuckets)+1
       ! Loop over all tetrahedrons that contain this ikibz as first element
       do jtetra=1,tetra_hash_count(ihash)
         ! if tetrahedron already exists add multiplicity
         if (tetra%unique_tetra(ihash)%indexes(1,jtetra)/=tetra_ibz(1)) cycle
         if (tetra%unique_tetra(ihash)%indexes(2,jtetra)/=tetra_ibz(2)) cycle
         if (tetra%unique_tetra(ihash)%indexes(3,jtetra)/=tetra_ibz(3)) cycle
         if (tetra%unique_tetra(ihash)%indexes(4,jtetra)/=tetra_ibz(4)) cycle
         tetra%unique_tetra(ihash)%indexes(0,jtetra) = tetra%unique_tetra(ihash)%indexes(0,jtetra)+1
         cycle tetra_loop1
       end do
       ! Otherwise store new tetrahedron
       tetra_hash_count(ihash) = tetra_hash_count(ihash)+1
       max_ntetra = size(tetra%unique_tetra(ihash)%indexes,2)
       ! The contents don't fit the array so I have to resize it
       if (tetra_hash_count(ihash)>max_ntetra) then
         ABI_MALLOC(indexes,(0:4,max_ntetra+TETRA_STEP))
         indexes(0:4,:max_ntetra) = tetra%unique_tetra(ihash)%indexes
         indexes(:,max_ntetra+1:) = 0
         ABI_MOVE_ALLOC(indexes,tetra%unique_tetra(ihash)%indexes)
       end if
       tetra%unique_tetra(ihash)%indexes(1:,tetra_hash_count(ihash)) = tetra_ibz(:)
       tetra%unique_tetra(ihash)%indexes(0, tetra_hash_count(ihash)) = 1
     end do tetra_loop1
   end do
   ! HM: The multiplicity of the tetrahedrons computed so far is wrong because we are using IBZ
   ! I compute the k-point multiplicity so I can fix this later on
   ! Only needed in blochl_weights* interface for looping over tetrahedra.
   ! in the onewk routines the weight is known outside
   ABI_CALLOC(tetra%ibz_multiplicity,(tetra%nkibz))
   do ikbz=1,nkpt_fullbz
     ikibz = bz2ibz(ikbz)
     tetra%ibz_multiplicity(ikibz) = tetra%ibz_multiplicity(ikibz) + 1
   end do
 case(2)
   min_idiag = 1
   ! For each k-point in the BZ generate the 6 tetrahedra that tesselate a microzone
   do ikbz=1,tetra%nkbz
     k1 = kpt_fullbz(:,ikbz)
     tetra_loop2: do itetra=1,6
       ! Determine tetrahedron
       do isummit=1,4
         ! Find the index of the neighbouring k-points in the BZ
         k2 = k1 + tetra_shifts_6(1,isummit,itetra,min_idiag)*klatt(:,1) + &
                   tetra_shifts_6(2,isummit,itetra,min_idiag)*klatt(:,2) + &
                   tetra_shifts_6(3,isummit,itetra,min_idiag)*klatt(:,3)
         ! Find full kpoint which is summit isummit of tetrahedron itetra around full kpt ikpt_full !
         !ikpt2 = octree_find(oct,k2,dist)
         !ikpt2 = octree_find_nearest_pbc(oct,k2,dist,shift)
         !if (dist>tol12) call exit(1)
         ikpt2 = kptrank_index(kptrank,k2)
         ! Find the index of those points in the BZ and IBZ
         tetra_ibz(isummit) = bz2ibz(ikpt2)
       end do
       ! Sort index of irr k-point edges (need this so the comparison works)
       call sort_4tetra_int(tetra_ibz)

       ! Store only unique tetrahedra
       ! Compute a very simple hash for each tetrahedron
       ihash = compute_hash(tetra,tetra_ibz) !mod(sum(tetra_ibz),tetra%nbuckets)+1
       ! Loop over all tetrahedrons that contain this ikibz as first element
       do jtetra=1,tetra_hash_count(ihash)
         ! if tetrahedron already exists add multiplicity
         if (tetra%unique_tetra(ihash)%indexes(1,jtetra)/=tetra_ibz(1)) cycle
         if (tetra%unique_tetra(ihash)%indexes(2,jtetra)/=tetra_ibz(2)) cycle
         if (tetra%unique_tetra(ihash)%indexes(3,jtetra)/=tetra_ibz(3)) cycle
         if (tetra%unique_tetra(ihash)%indexes(4,jtetra)/=tetra_ibz(4)) cycle
         tetra%unique_tetra(ihash)%indexes(0,jtetra) = tetra%unique_tetra(ihash)%indexes(0,jtetra)+1
         cycle tetra_loop2
       end do
       ! Otherwise store new tetrahedron
       tetra_hash_count(ihash) = tetra_hash_count(ihash)+1
       max_ntetra = size(tetra%unique_tetra(ihash)%indexes,2)
       ! The contents don't fit the array so I have to resize it
       if (tetra_hash_count(ihash)>max_ntetra) then
         ABI_MALLOC(indexes,(0:4,max_ntetra+TETRA_STEP))
         indexes(0:4,:max_ntetra) = tetra%unique_tetra(ihash)%indexes
         indexes(:,max_ntetra+1:) = 0
         ABI_MOVE_ALLOC(indexes,tetra%unique_tetra(ihash)%indexes)
       end if
       tetra%unique_tetra(ihash)%indexes(1:,tetra_hash_count(ihash)) = tetra_ibz(:)
       tetra%unique_tetra(ihash)%indexes(0, tetra_hash_count(ihash)) = 1
     end do tetra_loop2
   end do
 case default
   ierr = 1
   write(errorstring,*) 'Invalid option for the generation of tetrahedra,',ch10,&
                        'possible options are:',ch10,&
                        '1. Generate 24 tetrahedra per k-point',ch10,&
                        '2. Generate tetrahedra in the FBZ a map to IBZ (default)'
   return
 end select
 !ierr = octree_free(oct)
 ABI_FREE(tetra_hash_count)
 call destroy_kptrank(kptrank)

 ! Do some maintenance: free unused memory and count unique tetrahedra per IBZ point
 tetra%nunique_tetra = 0
 do ihash=1,tetra%nbuckets
   ! Count tetrahedra in this bucket
   ntetra = count(tetra%unique_tetra(ihash)%indexes(0,:)>0)
   tetra%nunique_tetra = tetra%nunique_tetra + ntetra
   ! Allocate array with right size
   ABI_MALLOC(indexes,(0:4,ntetra))
   indexes = tetra%unique_tetra(ihash)%indexes(:,:ntetra)
   ABI_MOVE_ALLOC(indexes,tetra%unique_tetra(ihash)%indexes)
 end do

 ! Sum the multiplicity
 ABI_MALLOC(tetra%tetra_count,(tetra%nkibz))
 ABI_MALLOC(tetra%tetra_total,(tetra%nkibz))
 tetra%tetra_count = 0
 tetra%tetra_total = 0
 do ihash=1,tetra%nbuckets
   ntetra = size(tetra%unique_tetra(ihash)%indexes,2)
   do itetra=1,ntetra
     tetra_mibz = tetra%unique_tetra(ihash)%indexes(:,itetra)
     do isummit=1,4
       ikibz = tetra_mibz(isummit)
       tetra%tetra_total(ikibz) = tetra%tetra_total(ikibz) + tetra_mibz(0)
       tetra%tetra_count(ikibz) = tetra%tetra_count(ikibz) + 1
     end do
   end do
 end do
 tetra%nibz_tetra = sum(tetra%tetra_count)

 ! HM: This was being allocated here, however this is only used when we loop over kpoints
 ! I will only allocate this memory if the htetra_get_onewk_* routines are called (lazy evaluation)
 !call htetra_init_mapping_ibz(tetra)

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

 contains
 integer function compute_hash(tetra,t) result(ihash)
   type(t_htetrahedron),intent(in) :: tetra
   integer,intent(in) :: t(4)
   ihash = mod(sum(t),tetra%nbuckets)+1
   ! TODO: should use a more general hash function that supports more buckets
   ! Something like:
   ! id = t(1)*nk3+t(2)*nk2+t(3)*nk1+t(4)
   ! where nk is the number of points in the IBZ
   ! Computing this leads to overflow so should use
   ! mod comutation operations
   ! (A + B) mod C = (A mod C + B mod C) mod C
   ! (A * B) mod C = (A mod C * B mod C) mod C
   ! A^B mod C = ( (A mod C)^B ) mod C
 end function compute_hash

end subroutine htetra_init
!!***

!----------------------------------------------------------------------

!!****f* m_htetrahedron/htetra_init_mapping_ibz
!! NAME
!! htetra_init_mapping_ibz
!!
!! FUNCTION
!!  The mapping to th IBZ is has its own allocation routine.
!!  I will only allocate this memory if the htetra_get_onewk_* routines are called (lazy evaluation)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine htetra_init_mapping_ibz(tetra)
 type(t_htetrahedron),intent(inout) :: tetra
 integer :: ikibz, itetra, isummit, ihash, ntetra
 integer :: tetra_count(tetra%nkibz),tetra_mibz(0:4)

 ! Only execute the following if not yet alocated
 if (allocated(tetra%ibz)) return

 ! Allocate IBZ to tetrahedron mapping
 ABI_MALLOC(tetra%ibz,(tetra%nkibz))
 do ikibz=1,tetra%nkibz
   ABI_MALLOC(tetra%ibz(ikibz)%indexes,(2,tetra%tetra_count(ikibz)))
 end do

 ! Create mapping from IBZ to unique tetrahedra
 tetra_count = 0
 do ihash=1,tetra%nbuckets
   ntetra = size(tetra%unique_tetra(ihash)%indexes,2)
   do itetra=1,ntetra
     tetra_mibz = tetra%unique_tetra(ihash)%indexes(:,itetra)
     do isummit=1,4
       ikibz = tetra_mibz(isummit)
       tetra_count(ikibz) = tetra_count(ikibz) + 1
       tetra%ibz(ikibz)%indexes(1,tetra_count(ikibz)) = ihash
       tetra%ibz(ikibz)%indexes(2,tetra_count(ikibz)) = itetra
     end do
   end do
 end do

end subroutine htetra_init_mapping_ibz
!!***


!----------------------------------------------------------------------

!!****f* m_htetrahedron/htetra_get_ibz
!! NAME
!! htetra_get_ibz
!!
!! FUNCTION
!!  Get the itetra tetrahedron contributing to the ikibz k-point
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

pure subroutine htetra_get_ibz(tetra,ikibz,itetra,tetra_mibz)
 type(t_htetrahedron), intent(in) :: tetra
 integer,intent(in) :: ikibz, itetra
 integer,intent(out) :: tetra_mibz(0:4)
 integer :: ihash, jtetra

 ihash  = tetra%ibz(ikibz)%indexes(1,itetra)
 jtetra = tetra%ibz(ikibz)%indexes(2,itetra)
 tetra_mibz = tetra%unique_tetra(ihash)%indexes(:,jtetra)
end subroutine htetra_get_ibz
!!***

!----------------------------------------------------------------------

!!****f* m_htetrahedron/htetra_print
!! NAME
!! htetra_print
!!
!! FUNCTION
!!  write information about the tetrahedra object
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine htetra_print(self)

 type(t_htetrahedron), intent(in) :: self
 real(dp) :: total_size, unique_tetra_size, ibz_pointer_size

 unique_tetra_size = self%nunique_tetra*5*four/1024/1024
 total_size        = unique_tetra_size
 write(std_out,'(a,i12)')     'unique_tetra      ', self%nunique_tetra
 write(std_out,'(a,f12.1,a)') 'unique_tetra_size ', unique_tetra_size, ' [Mb]'
 if (allocated(self%ibz)) then
   ibz_pointer_size  = self%nibz_tetra*2*four/1024/1024
   write(std_out,'(a,i12)')     'ibz_tetra         ', self%nibz_tetra
   write(std_out,'(a,f12.1,a)') 'ibz_tetra_size    ', ibz_pointer_size, ' [Mb]'
   total_size = total_size + ibz_pointer_size
 end if
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
 integer :: ikibz,ihash

 ABI_SFREE(tetra%tetra_count)
 ABI_SFREE(tetra%tetra_total)
 ABI_SFREE(tetra%ibz_multiplicity)

 if (allocated(tetra%unique_tetra)) then
   do ihash=1,tetra%nbuckets
     ABI_SFREE(tetra%unique_tetra(ihash)%indexes)
   end do
   ABI_FREE(tetra%unique_tetra)
 end if

 if (allocated(tetra%ibz)) then
   do ikibz=1,tetra%nkibz
     ABI_SFREE(tetra%ibz(ikibz)%indexes)
   end do
   ABI_FREE(tetra%ibz)
 end if

end subroutine htetra_free
!!***

!----------------------------------------------------------------------

!!****f* m_htetrahedron/get_onetetra_blochl
!! NAME
!! get_onetetra_blochl
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

pure subroutine get_onetetra_blochl(eig,energies,nene,bcorr,tweight,dweight)

!Arguments ------------------------------------
!scalars
 integer,intent(in)   :: nene,bcorr
 real(dp) ,intent(in) :: energies(nene)
!arrays
 real(dp), intent(out) :: tweight(4, nene)
 real(dp), intent(out) :: dweight(4, nene)
 real(dp),intent(in)  :: eig(4)

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
 real(dp) :: bcorr_fact

! *********************************************************************

 ! This is output
 tweight = zero; dweight = zero

 ! all notations are from Blochl PRB 49 16223 [[cite:Bloechl1994a]] Appendix B
 e1 = eig(1)
 e2 = eig(2)
 e3 = eig(3)
 e4 = eig(4)
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
     cc = inv_e21*inv_e31*inv_e41*deleps1**3
     tweight(1,ieps) = cc*(4.d0-deleps1*invepsum)
     tweight(2,ieps) = cc*deleps1*inv_e21
     tweight(3,ieps) = cc*deleps1*inv_e31
     tweight(4,ieps) = cc*deleps1*inv_e41

     ! Delta
     dccde_pre = 3.d0*inv_e21*inv_e31*inv_e41
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
     cc1_pre = inv_e31*inv_e41
     cc2_pre = inv_e41*inv_e32*inv_e31
     cc3_pre = inv_e42*inv_e32*inv_e41

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
     dcc1de = cc1_pre*(2.d0*deleps1)
     dcc2de = cc2_pre*(    -deleps1*deleps2+deleps1*deleps3+deleps2*deleps3)
     dcc3de = cc3_pre*(2.d0*deleps2*deleps4-deleps2*deleps2)
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
     cc_pre = inv_e41*inv_e42*inv_e43
     cc = cc_pre*deleps4**3
     tweight(1,ieps) = one - deleps4*cc*inv_e41
     tweight(2,ieps) = one - deleps4*cc*inv_e42
     tweight(3,ieps) = one - deleps4*cc*inv_e43
     tweight(4,ieps) = one - cc*(4.d0-deleps4*invepsum)

     ! Delta
     dccde = -3.d0*cc_pre*deleps4**2
     dccde_tmp = dccde*deleps4 - cc
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
     tweight(:,ieps:) = one

     ! Delta unchanged by this tetrahedron
     exit
   end if

   !
   !  if we have a fully degenerate tetrahedron,
   !  1) the tweight is a Heaviside (step) function, which is correct above, but
   !  2) the dweight should contain a Dirac function
   !
 end do

end subroutine get_onetetra_blochl
!!***

!!****f* m_htetrahedron/get_ontetra_lambinvigneron
!! NAME
!! get_ontetra_lambinvigneron
!!
!! FUNCTION
!!  Compute the complex weights according to:
!!  P. Lambin and J.P. Vigneron, Phys. Rev. B 29, 3430 (1984)
!!  This routine is adapted from tdep where it was implemented
!!  by Olle Hellman, all credits go to him
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

pure subroutine get_ontetra_lambinvigneron(eig,z,cw)
    ! dispersion values at the corners of the tetrahedron
    real(dp), intent(in) :: eig(4)
    ! energy to evaulate the weights at
    complex(dp), intent(in) :: z
    ! complex weights
    complex(dp), intent(out) :: cw(4)
    complex(dp) :: EZ1,EZ2,EZ3,EZ4
    real(dp) :: Emin,Emax,Zdist
    real(dp) :: E12,E13,E14,E23,E24,E34
    real(dp) :: a,b,c,d,e,f

    cw = zero

    ! Min and max energy
    Emin=eig(1)
    Emax=eig(4)

    ! First the complex energy differences
    EZ1=z-eig(1)
    EZ2=z-eig(2)
    EZ3=z-eig(3)
    EZ4=z-eig(4)
    ! Smallest distance |z-Ei|, to determine wether I should switch to the
    ! asymptotic behavior for numerical stability.
    Zdist=huge(Zdist)
    Zdist=min(Zdist,abs(EZ2))
    Zdist=min(Zdist,abs(EZ3))
    Zdist=min(Zdist,abs(EZ4))
    !@TODO add asymptotic thing with continued fractions

    ! Then the energy differences, for the coefficients. Must always be positive, I hope.
    E12=eig(2)-eig(1)
    E13=eig(3)-eig(1)
    E14=eig(4)-eig(1)
    E23=eig(3)-eig(2)
    E24=eig(4)-eig(2)
    E34=eig(4)-eig(3)
    a=zero; if ( E12 .gt. tol14 ) a=one/E12
    b=zero; if ( E13 .gt. tol14 ) b=one/E13
    c=zero; if ( E14 .gt. tol14 ) c=one/E14
    d=zero; if ( E23 .gt. tol14 ) d=one/E23
    e=zero; if ( E24 .gt. tol14 ) e=one/E24
    f=zero; if ( E34 .gt. tol14 ) f=one/E34

    ! Now get the actual weights
    ! e1=e2=e3=e4
    if ( E12+E23+E34 .lt. tol14 ) then
        cw(1)=0.25_dp/EZ1
        cw(2)=0.25_dp/EZ2
        cw(3)=0.25_dp/EZ3
        cw(4)=0.25_dp/EZ4
    !    e2=e3=e4
    elseif ( E23+E34 .lt. tol14 ) then
        cw(1)=-a - (3*a**2*EZ2)*half + 3*a**3*EZ1*EZ2 + 3*a**4*EZ1**2*EZ2*Log(EZ2/EZ1)
        cw(2)=-a - (3*a**2*EZ2)*half + 3*a**3*EZ1*EZ2 + 3*a**4*EZ1**2*EZ2*Log(EZ2/EZ1)
        cw(3)=-a - (3*a**2*EZ2)*half + 3*a**3*EZ1*EZ2 + 3*a**4*EZ1**2*EZ2*Log(EZ2/EZ1)
        cw(4)=-a*third + (a**2*EZ1)*half - a**3*EZ1**2 + a**4*EZ1**3*Log(EZ2/EZ1)
    ! e1=e2=e3
    elseif ( E12+E23 .lt. tol14 ) then
        cw(1)=f*third - (EZ4*f**2)*half + EZ4**2*f**3 + EZ4**3*f**4*Log(EZ4/EZ3)
        cw(2)=f*third - (EZ4*f**2)*half + EZ4**2*f**3 + EZ4**3*f**4*Log(EZ4/EZ3)
        cw(3)=f*third - (EZ4*f**2)*half + EZ4**2*f**3 + EZ4**3*f**4*Log(EZ4/EZ3)
        cw(4)=-f + (3*EZ3*f**2)*half - 3*EZ3*EZ4*f**3 + 3*EZ3*EZ4**2*f**4*Log(EZ3/EZ4)
    ! e1=e2 e3=e4
    elseif ( E12+E34 .lt. tol14 ) then
        cw(1)=-d - (3*d**2*EZ2)*half + 3*d**3*EZ2*EZ3 + 3*d**4*EZ2*EZ3**2*Log(EZ2/EZ3)
        cw(2)=-d - (3*d**2*EZ2)*half + 3*d**3*EZ2*EZ3 + 3*d**4*EZ2*EZ3**2*Log(EZ2/EZ3)
        cw(3)=d - (3*d**2*EZ3)*half - 3*d**3*EZ2*EZ3 + 3*d**4*EZ2**2*EZ3*Log(EZ3/EZ2)
        cw(4)=d - (3*d**2*EZ3)*half - 3*d**3*EZ2*EZ3 + 3*d**4*EZ2**2*EZ3*Log(EZ3/EZ2)
    !       e3=e4
    elseif ( E34 .lt. tol14 ) then
        cw(1)=-(a*b**2*EZ1**2*(-1 + (a*EZ2 + 2*b*EZ3)*Log(EZ1))) + a**2*d**2*EZ2**3*Log(EZ2) - &
                b**2*d*EZ3**2*(1 + (2*b*EZ1 + d*EZ2)*Log(EZ3))
        cw(2)=a**2*b**2*EZ1**3*Log(EZ1) - a*d**2*EZ2**2*(1 + (a*EZ1 - 2*d*EZ3)*Log(EZ2)) - &
              b*d**2*EZ3**2*(1 + (b*EZ1 + 2*d*EZ2)*Log(EZ3))
        cw(3)=a*b**3*EZ1**3*Log(EZ1) - a*d**3*EZ2**3*Log(EZ2) + b*d*EZ3*(half + b*EZ1 + d*EZ2 + &
             (b**2*EZ1**2 + b*d*EZ1*EZ2 + d**2*EZ2**2)*Log(EZ3))
        cw(4)=a*b**3*EZ1**3*Log(EZ1) - a*d**3*EZ2**3*Log(EZ2) + b*d*EZ3*(half + b*EZ1 + d*EZ2 + &
             (b**2*EZ1**2 + b*d*EZ1*EZ2 + d**2*EZ2**2)*Log(EZ3))
    !    e2=e3
    elseif ( E23 .lt. tol14 ) then
        cw(1)=-(a**2*c*EZ1**2*(-1 + (2*a*EZ2 + c*EZ4)*Log(EZ1))) + &
                a**2*e*EZ2**2*(1 + (2*a*EZ1 - e*EZ4)*Log(EZ2)) + c**2*e**2*EZ4**3*Log(EZ4)
        cw(2)=a**3*c*EZ1**3*Log(EZ1) - &
              a*e*EZ2*(half + a*EZ1 - e*EZ4 + (a**2*EZ1**2 - a*e*EZ1*EZ4 + e**2*EZ4**2)*Log(EZ2)) + c*e**3*EZ4**3*Log(EZ4)
        cw(3)=a**3*c*EZ1**3*Log(EZ1) - &
              a*e*EZ2*(half + a*EZ1 - e*EZ4 + (a**2*EZ1**2 - a*e*EZ1*EZ4 + e**2*EZ4**2)*Log(EZ2)) + c*e**3*EZ4**3*Log(EZ4)
        cw(4)=a**2*c**2*EZ1**3*Log(EZ1) - &
              a*e**2*EZ2**2*(1 + (a*EZ1 - 2*e*EZ4)*Log(EZ2)) - c*e**2*EZ4**2*(1 + (c*EZ1 + 2*e*EZ2)*Log(EZ4))
    ! e1=e2
    elseif ( E12 .lt. tol14 ) then
        cw(1)=b*c*EZ1*(half - b*EZ3 - c*EZ4 + (b**2*EZ3**2 + b*c*EZ3*EZ4 + c**2*EZ4**2)*Log(EZ1)) - &
              b**3*EZ3**3*f*Log(EZ3) + c**3*EZ4**3*f*Log(EZ4)
        cw(2)=b*c*EZ1*(half - b*EZ3 - c*EZ4 + (b**2*EZ3**2 + b*c*EZ3*EZ4 + c**2*EZ4**2)*Log(EZ1)) - &
              b**3*EZ3**3*f*Log(EZ3) + c**3*EZ4**3*f*Log(EZ4)
        cw(3)=-(b**2*c*EZ1**2*(-1 + (2*b*EZ3 + c*EZ4)*Log(EZ1))) + &
                b**2*EZ3**2*f*(1 + (2*b*EZ1 - EZ4*f)*Log(EZ3)) + c**2*EZ4**3*f**2*Log(EZ4)
        cw(4)=-(b*c**2*EZ1**2*(-1 + (b*EZ3 + 2*c*EZ4)*Log(EZ1))) + &
                b**2*EZ3**3*f**2*Log(EZ3) - c**2*EZ4**2*f*(1 + (2*c*EZ1 + EZ3*f)*Log(EZ4))
    ! e1<e2<e3<e4
    else
        cw(1)=a**2*d*e*EZ2**3*Log(EZ2/EZ1) - b**2*d*EZ3**3*f*Log(EZ3/EZ1) + c*(a*b*EZ1**2 + c*e*EZ4**3*f*Log(EZ4/EZ1))
        cw(2)=a**2*b*c*EZ1**3*Log(EZ1/EZ2) - b*d**2*EZ3**3*f*Log(EZ3/EZ2) + e*(-(a*d*EZ2**2) + c*e*EZ4**3*f*Log(EZ4/EZ2))
        cw(3)=a*b**2*c*EZ1**3*Log(EZ1/EZ3) - a*d**2*e*EZ2**3*Log(EZ2/EZ3) + f*(b*d*EZ3**2 + c*e*EZ4**3*f*Log(EZ4/EZ3))
        cw(4)=a*b*c**2*EZ1**3*Log(EZ1/EZ4) - a*d*e**2*EZ2**3*Log(EZ2/EZ4) + f*(-(c*e*EZ4**2) + b*d*EZ3**3*f*Log(EZ3/EZ4))
    endif

    ! HM:check this
    cw = cw * two

end subroutine get_ontetra_lambinvigneron
!!***

!!****f* m_htetrahedron/get_ontetratra_lambinvigneron_imag
!! NAME
!! get_ontetratra_lambinvigneron_imag
!!
!! FUNCTION
!!  Compute the complex weights according to:
!!  P. Lambin and J.P. Vigneron, Phys. Rev. B 29, 3430 (1984)
!!  This routine is adapted from tdep where it was implemented
!!  by Olle Hellman, all credits go to him
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

pure subroutine get_ontetetra_lambinvigneron_imag(eig,energies,nene,wt)
 ! dispersion values at the corners of the tetrahedron
 real(dp), intent(in), dimension(4) :: eig
 ! number of energies
 integer, intent(in) :: nene
 ! energy to evaulate the weights at
 real(dp), intent(in) :: energies(nene)
 ! integration weights
 real(dp), intent(out) :: wt(4,nene)

 integer :: ie
 real(dp) :: z
 real(dp) :: EZ1,EZ2,EZ3,EZ4
 real(dp) :: Emin,Emax
 real(dp) :: E12,E13,E14,E23,E24,E34
 real(dp) :: a,b,c,d,e,f,ff0,ff1,ff2,ff3,gg0,gg1,gg2,gg3,hh0,hh1,hh2,hh3,ii0,ii1,ii2,ii3

 wt = zero
 Emin = eig(1)
 Emax = eig(4)

 do ie=1,nene
   z = energies(ie)
   if (z<eig(1)) then ! e<e1<e2<e3<e4
     cycle
   else if (z .lt. eig(2)) then ! e1<e<e2<e3<e4
     EZ1=z-eig(1)
     EZ2=z-eig(2)
     EZ3=z-eig(3)
     EZ4=z-eig(4)
     E12=eig(2)-eig(1)
     E13=eig(3)-eig(1)
     E14=eig(4)-eig(1)
     a=one/E12
     b=one/E13
     c=one/E14
     wt(1,ie)=a*b*c*EZ1**2*(-a*EZ2 - b*EZ3 - c*EZ4)
     wt(2,ie)=a**2*b*c*EZ1**3
     wt(3,ie)=a*b**2*c*EZ1**3
     wt(4,ie)=a*b*c**2*EZ1**3
     cycle
   else if (z .lt. eig(3)) then ! e1<e2<e<e3<e4
     EZ1=z-eig(1)
     EZ2=z-eig(2)
     EZ3=z-eig(3)
     EZ4=z-eig(4)
     E13=eig(3)-eig(1)
     E14=eig(4)-eig(1)
     E23=eig(3)-eig(2)
     E24=eig(4)-eig(2)
     b=one/E13
     c=one/E14
     d=one/E23
     e=one/E24
     ff0=-b**2*EZ3
     ff2=-c**2*EZ4
     gg0=-d**2*EZ3
     gg2=-e**2*EZ4
     hh0=d**2*EZ2
     hh2=b**2*EZ1
     ii0=e**2*EZ2
     ii2=c**2*EZ1
     ff1=-c*d*EZ1*EZ3-d*e*EZ2*EZ3-c*e*EZ1*EZ4
     ff3=-b*d*EZ1*EZ3-b*e*EZ1*EZ4-d*e*EZ2*EZ4
     gg1=-b*c*EZ1*EZ3-b*e*EZ2*EZ3-c*e*EZ2*EZ4
     gg3=-b*d*EZ2*EZ3-b*c*EZ1*EZ4-c*d*EZ2*EZ4
     hh1=-b*c*EZ1*EZ3-b*e*EZ2*EZ3-c*e*EZ2*EZ4
     hh3=-c*d*EZ1*EZ3-d*e*EZ2*EZ3-c*e*EZ1*EZ4
     ii1=-b*d*EZ2*EZ3-b*c*EZ1*EZ4-c*d*EZ2*EZ4
     ii3=-b*d*EZ1*EZ3-b*e*EZ1*EZ4-d*e*EZ2*EZ4
     wt(1,ie)=half*(ff0*ff1+ff2*ff3)
     wt(2,ie)=half*(gg0*gg1+gg2*gg3)
     wt(3,ie)=half*(hh0*hh1+hh2*hh3)
     wt(4,ie)=half*(ii0*ii1+ii2*ii3)
     cycle
   else if (z .lt. eig(4)) then ! e1<e2<e3<e<e4
     EZ1=z-eig(1)
     EZ2=z-eig(2)
     EZ3=z-eig(3)
     EZ4=z-eig(4)
     E14=eig(4)-eig(1)
     E24=eig(4)-eig(2)
     E34=eig(4)-eig(3)
     c=one/E14
     e=one/E24
     f=one/E34
     wt(1,ie)=-(c**2*e*EZ4**3*f)
     wt(2,ie)=-(c*e**2*EZ4**3*f)
     wt(3,ie)=-(c*e*EZ4**3*f**2)
     wt(4,ie)=c*e*EZ4**2*f*(c*EZ1 + e*EZ2 + EZ3*f)
     cycle
   else
     exit
   end if
 end do
 wt = wt*4.0_dp

end subroutine get_ontetetra_lambinvigneron_imag
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

subroutine htetra_get_onewk_wvals(tetra, ik_ibz, opt, nw, wvals, max_occ, nkibz, eig_ibz, weights)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_ibz,nw,nkibz,opt
 real(dp) ,intent(in) :: max_occ
 type(t_htetrahedron), intent(inout) :: tetra
!arrays
 real(dp),intent(in) :: wvals(nw)
 real(dp),intent(in) :: eig_ibz(nkibz)
 real(dp),intent(out) :: weights(nw, 2)

!Local variables-------------------------------
!scalars
 integer  :: itetra,isummit,tetra_count,tetra_total
 real(dp) :: tweight
!arrays
 integer  :: ind_ibz(4),tetra_mibz(0:4)
 real(dp) :: eig(4)
 real(dp) :: tweight_tmp(4,nw),dweight_tmp(4,nw)

! *********************************************************************

 weights = zero
 ! lazy evaluation of the mapping from k-points to tetrahedra
 if (.not.allocated(tetra%ibz)) call htetra_init_mapping_ibz(tetra)

 ! For each tetrahedron that belongs to this k-point
 tetra_count = tetra%tetra_count(ik_ibz)
 tetra_total = tetra%tetra_total(ik_ibz)
 do itetra=1,tetra_count

   call htetra_get_ibz(tetra,ik_ibz,itetra,tetra_mibz)
   tweight = one*tetra_mibz(0)/tetra_total
   do isummit=1,4
     ! Get mapping of each summit to eig_ibz
     ind_ibz(isummit) = tetra_mibz(isummit)
     eig(isummit) = eig_ibz(ind_ibz(isummit))
   end do

   ! Sort energies before calling get_onetetra_blochl
   !call sort_tetra(4, eig, ind_ibz, tol14)
   call sort_4tetra(eig, ind_ibz)

   ! HM: Here we should only compute what we will use!
   select case (opt)
   case(0:1)
     call get_onetetra_blochl(eig, wvals, nw, opt, tweight_tmp, dweight_tmp)
   case(2)
     call get_ontetetra_lambinvigneron_imag(eig, wvals, nw, dweight_tmp)
     tweight_tmp = zero
   end select

   ! Accumulate contributions to ik_ibz (there might be multiple vertexes that map onto ik_ibz)
   do isummit=1,4
     if (ind_ibz(isummit) /= ik_ibz) cycle
     weights(:,1) = weights(:,1) + dweight_tmp(isummit,:)*tweight*max_occ
     weights(:,2) = weights(:,2) + tweight_tmp(isummit,:)*tweight*max_occ
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
 type(t_htetrahedron), intent(inout) :: tetra
 real(dp) ,intent(in) :: enemin,enemax,max_occ
!arrays
 real(dp),intent(in) :: eig_ibz(nkibz)
 real(dp),intent(out) :: weights(nw,2)

!Local variables-------------------------------
!scalars
 real(dp) :: wvals(nw)

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

subroutine htetra_get_onewk_wvals_zinv(tetra, ik_ibz, nz, zvals, max_occ, nkibz, eig_ibz, opt, cweights)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_ibz,nz,nkibz,opt
 real(dp) ,intent(in) :: max_occ
 type(t_htetrahedron), intent(inout) :: tetra
!arrays
 complex(dp),intent(in) :: zvals(nz)
 real(dp),intent(in) :: eig_ibz(nkibz)
 complex(dp),intent(out) :: cweights(nz)

!Local variables-------------------------------
!scalars
 integer  :: itetra,isummit,tetra_total,tetra_count,iz
 real(dp) :: tweight
!arrays
 integer  :: ind_ibz(4),tetra_mibz(0:4)
 real(dp) :: eig(4)
 complex(dp) :: verm(4), cw(4), verli(4)
! *********************************************************************

 cweights = zero
 ! lazy evaluation of the mapping from k-points to tetrahedra
 if (.not.allocated(tetra%ibz)) call htetra_init_mapping_ibz(tetra)

 ! For each tetrahedron that belongs to this k-point
 tetra_count = tetra%tetra_count(ik_ibz)
 tetra_total = tetra%tetra_total(ik_ibz)
 do itetra=1,tetra_count

   call htetra_get_ibz(tetra,ik_ibz,itetra,tetra_mibz)
   tweight = one*tetra_mibz(0)/tetra_total
   do isummit=1,4
     ! Get mapping of each summit to eig_ibz
     ind_ibz(isummit) = tetra_mibz(isummit)
     eig(isummit) = eig_ibz(ind_ibz(isummit))
   end do

   ! Loop over frequencies
   do iz=1,nz
     select case(opt)
     case(1)
       verm = zvals(iz) - eig
       call SIM0TWOI(cw, VERLI, VERM)
     case(2)
       call get_ontetra_lambinvigneron(eig,zvals(iz),cw)
     end select

     do isummit=1,4
       if (ind_ibz(isummit) /= ik_ibz) cycle
       cweights(iz) = cweights(iz) + cw(isummit) *tweight*max_occ
       ! HM: This exit is important, avoids summing the same contribution more than once
       exit
     end do
   end do
 end do ! itetra

end subroutine htetra_get_onewk_wvals_zinv
!!***

!----------------------------------------------------------------------

!!****f* m_htetrahedron/htetra_get_delta_mask
!! NAME
!!  htetra_get_delta_mask
!!
!! FUNCTION
!!  Get a mask for the kpoints where the delta is finite
!!

subroutine htetra_get_delta_mask(tetra,eig_ibz,wvals,nw,nkpt,kmask,comm)
!Arguments
 integer,intent(in) :: nw,nkpt,comm
 type(t_htetrahedron), intent(in) :: tetra
 real(dp),intent(in) :: wvals(nw)
 real(dp),intent(in) :: eig_ibz(nkpt)
 integer,intent(out) :: kmask(nkpt)

!Local variables-------------------------------
 integer :: ik_ibz,nprocs,my_rank,ierr
 integer :: tetra_count, itetra, isummit, ihash
 integer :: contrib
 real(dp) :: emin,emax
 integer :: ind_ibz(4)
 real(dp) :: eig(4)

 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 kmask = 0
 ! For each bucket of tetrahedra
 do ihash=1,tetra%nbuckets
   if (mod(ihash,nprocs) /= my_rank) cycle

   ! For each tetrahedron
   tetra_count = size(tetra%unique_tetra(ihash)%indexes,2)
   do itetra=1,tetra_count

     ! Get mapping of each summit to eig_ibz
     do isummit=1,4
       ind_ibz(isummit) = tetra%unique_tetra(ihash)%indexes(isummit,itetra)
       eig(isummit) = eig_ibz(ind_ibz(isummit))
     end do

     ! Determine the energy range of the tetrahedra
     emin = minval(eig)
     emax = maxval(eig)

     ! Check if any value in wvals is betwen emin and emax
     contrib = 0; if (any(emin<wvals.and.wvals<emax)) contrib = 1

     ! Compute the union
     do isummit=1,4
       ik_ibz = ind_ibz(isummit)
       kmask(ik_ibz) = kmask(ik_ibz) + contrib
     end do
   end do ! itetra
 end do

 call xmpi_sum(kmask, comm, ierr)

end subroutine htetra_get_delta_mask
!!***

!----------------------------------------------------------------------

!!****f* m_htetrahedron/htetra_wvals_weights
!! NAME
!!  htetra_wvals_weights
!!
!! FUNCTION
!!   Emulates the behaviour of the previous tetrahedron implementation but
!!   taking a list of energies as input.
!!   HM: I find that in many routines its better to change the implementation
!!   and accumulate the tetrahedron weights in the same way as the
!!   gaussian smearing weights using htetra_get_onewk_wvals. However this requires
!!   some refactoring of the code. I provide this routine to make it easier
!!   to transition to the new tetrahedron implementation without refactoring.
!!   Looping over tetrahedra (i.e. using tetra_blochl_weights) is currently faster
!!   than looping over k-points.
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

subroutine htetra_wvals_weights(tetra,eig_ibz,nw,wvals,max_occ,nkpt,opt,tweight,dweight,comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nw,nkpt,opt,comm
 type(t_htetrahedron), intent(in) :: tetra
 real(dp) ,intent(in) :: max_occ
!arrays
 real(dp),intent(in) :: eig_ibz(nkpt)
 real(dp),intent(out) :: dweight(nw,nkpt),tweight(nw,nkpt)

!Local variables-------------------------------
!scalars
 integer :: ik_ibz,multiplicity,nprocs,my_rank,ierr
 integer :: tetra_count, itetra, isummit, ihash
!arrays
 integer :: ind_ibz(4)
 real(dp) :: eig(4)
 real(dp) :: wvals(nw)
 real(dp) :: dweight_tmp(4,nw),tweight_tmp(4,nw)

! *********************************************************************

 tweight = zero; dweight = zero
 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 ! For each bucket of tetrahedra
 do ihash=1,tetra%nbuckets
   if (mod(ihash,nprocs) /= my_rank) cycle

   ! For each tetrahedron
   tetra_count = size(tetra%unique_tetra(ihash)%indexes,2)
   do itetra=1,tetra_count

     ! Get mapping of each summit to eig_ibz
     do isummit=1,4
       ind_ibz(isummit) = tetra%unique_tetra(ihash)%indexes(isummit,itetra)
       eig(isummit) = eig_ibz(ind_ibz(isummit))
     end do

     ! Sort energies before calling get_onetetra_blochl
     call sort_4tetra(eig, ind_ibz)

     ! Get tetrahedron weights
     select case (opt)
     case(0:1)
       call get_onetetra_blochl(eig, wvals, nw, opt, tweight_tmp, dweight_tmp)
     case(2)
       call get_ontetetra_lambinvigneron_imag(eig, wvals, nw, dweight_tmp)
       tweight_tmp = zero
     end select

     ! Acumulate the contributions
     multiplicity = tetra%unique_tetra(ihash)%indexes(0,itetra)
     do isummit=1,4
       ik_ibz = ind_ibz(isummit)
       dweight(:,ik_ibz) = dweight(:,ik_ibz) + dweight_tmp(isummit,:)*multiplicity*max_occ
       tweight(:,ik_ibz) = tweight(:,ik_ibz) + tweight_tmp(isummit,:)*multiplicity*max_occ
     end do
   end do ! itetra
 end do

 ! Rescale weights
 select case(tetra%opt)
 case(1)
   do ik_ibz=1,tetra%nkibz
     dweight(:,ik_ibz) = dweight(:,ik_ibz)*tetra%ibz_multiplicity(ik_ibz)/tetra%tetra_total(ik_ibz)/tetra%nkbz
     tweight(:,ik_ibz) = tweight(:,ik_ibz)*tetra%ibz_multiplicity(ik_ibz)/tetra%tetra_total(ik_ibz)/tetra%nkbz
   end do
 case(2)
   dweight = dweight*tetra%vv/4.0_dp
   tweight = tweight*tetra%vv/4.0_dp
 end select

 call xmpi_sum(dweight, comm, ierr)
 call xmpi_sum(tweight, comm, ierr)

end subroutine htetra_wvals_weights
!!***

!----------------------------------------------------------------------

!!****f* m_htetrahedron/htetra_wvals_weights_delta
!! NAME
!!  htetra_wvals_weights_delta
!!
!! FUNCTION
!!  Same as above but computing only delta for performance and memory
!!  HM: Should find a clean way to avoid copy paste routine
!!

subroutine htetra_wvals_weights_delta(tetra,eig_ibz,nw,wvals,max_occ,nkpt,opt,dweight,comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nw,nkpt,opt,comm
 type(t_htetrahedron), intent(in) :: tetra
 real(dp) ,intent(in) :: max_occ
!arrays
 real(dp),intent(in) :: eig_ibz(nkpt)
 real(dp),intent(out) :: dweight(nw,nkpt)

!Local variables-------------------------------
!scalars
 integer :: ik_ibz,multiplicity,nprocs,my_rank,ierr
 integer :: tetra_count, itetra, isummit, ihash
!arrays
 integer :: ind_ibz(4)
 real(dp) :: eig(4)
 real(dp) :: wvals(nw)
 real(dp) :: dweight_tmp(4,nw),tweight_tmp(4,nw)

! *********************************************************************

 dweight = zero
 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 ! For each bucket of tetrahedra
 do ihash=1,tetra%nbuckets
   if (mod(ihash,nprocs) /= my_rank) cycle

   ! For each tetrahedron
   tetra_count = size(tetra%unique_tetra(ihash)%indexes,2)
   do itetra=1,tetra_count

     ! Get mapping of each summit to eig_ibz
     do isummit=1,4
       ind_ibz(isummit) = tetra%unique_tetra(ihash)%indexes(isummit,itetra)
       eig(isummit) = eig_ibz(ind_ibz(isummit))
     end do

     ! Sort energies before calling get_onetetra_blochl
     call sort_4tetra(eig, ind_ibz)

     ! Get tetrahedron weights
     select case (opt)
     case(0:1)
       call get_onetetra_blochl(eig, wvals, nw, opt, tweight_tmp, dweight_tmp)
     case(2)
       call get_ontetetra_lambinvigneron_imag(eig, wvals, nw, dweight_tmp)
     end select

     ! Acumulate the contributions
     multiplicity = tetra%unique_tetra(ihash)%indexes(0,itetra)
     do isummit=1,4
       ik_ibz = ind_ibz(isummit)
       dweight(:,ik_ibz) = dweight(:,ik_ibz) + dweight_tmp(isummit,:)*multiplicity*max_occ
     end do
   end do ! itetra
 end do

 ! Rescale weights
 select case(tetra%opt)
 case(1)
   do ik_ibz=1,tetra%nkibz
     dweight(:,ik_ibz) = dweight(:,ik_ibz)*tetra%ibz_multiplicity(ik_ibz)/tetra%tetra_total(ik_ibz)/tetra%nkbz
   end do
 case(2)
   dweight = dweight*tetra%vv/4.0_dp
 end select

 call xmpi_sum(dweight, comm, ierr)

end subroutine htetra_wvals_weights_delta
!!***

!----------------------------------------------------------------------

!!****f* m_htetrahedron/htetra_blochl_weights
!! NAME
!!  htetra_blochl_weights
!!
!! FUNCTION
!!   Emulates the behaviour of the previous tetrahedron implementation.
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
 real(dp) :: wvals(nw)

! *********************************************************************

 wvals = linspace(enemin,enemax,nw)
 call htetra_wvals_weights(tetra,eig_ibz,nw,wvals,max_occ,nkpt,bcorr,tweight,dweight,comm)

end subroutine htetra_blochl_weights
!!***

!----------------------------------------------------------------------

!!****f* m_htetrahedron/htetra_blochl_weights_wvals_zinv
!! NAME
!!  htetra_blochl_weights_wvals_zinv
!!
!! FUNCTION
!!   The same as htetra_get_onewk_wvals_zinv but looping over tetrahedra
!!   which is more efficient
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

subroutine htetra_weights_wvals_zinv(tetra,eig_ibz,nz,zvals,max_occ,nkpt,opt,cweight,comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nz,nkpt,opt,comm
 type(t_htetrahedron), intent(in) :: tetra
 real(dp) ,intent(in) :: max_occ
!arrays
 real(dp) ,intent(in) :: eig_ibz(nkpt)
 complex(dp),intent(in)  :: zvals(nz)
 complex(dp),intent(out) :: cweight(nz,nkpt)

!Local variables-------------------------------
!scalars
 integer :: ik_ibz,iz,multiplicity,nprocs,my_rank,ierr
 integer :: tetra_count, itetra, isummit, ihash
!arrays
 integer :: ind_ibz(4)
 real(dp) :: eig(4)
 complex(dp) :: cw(4), verli(4), verm(4)
! *********************************************************************

 cweight = zero
 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 ! For each bucket of tetrahedra
 do ihash=1,tetra%nbuckets
   if (mod(ihash,nprocs) /= my_rank) cycle

   ! For each tetrahedron that belongs to this k-point
   tetra_count = size(tetra%unique_tetra(ihash)%indexes,2)
   do itetra=1,tetra_count

     ! Get mapping of each summit to eig_ibz
     do isummit=1,4
       ind_ibz(isummit) = tetra%unique_tetra(ihash)%indexes(isummit,itetra)
       eig(isummit) = eig_ibz(ind_ibz(isummit))
     end do

     ! Get multiplicity
     multiplicity = tetra%unique_tetra(ihash)%indexes(0,itetra)

     ! Loop over frequencies
     do iz=1,nz
       ! Get tetrahedron weights
       select case(opt)
       case(1)
         verm = zvals(iz) - eig
         call SIM0TWOI(cw, VERLI, VERM)
       case(2)
         call get_ontetra_lambinvigneron(eig,zvals(iz),cw)
       end select

       ! Acumulate the contributions
       do isummit=1,4
         ik_ibz = ind_ibz(isummit)
         cweight(iz,ik_ibz) = cweight(iz,ik_ibz) + cw(isummit)*multiplicity*max_occ
       end do
     end do ! iz
   end do ! itetra
 end do

 ! Rescale weights
 select case(tetra%opt)
 case(1)
   do ik_ibz=1,tetra%nkibz
     cweight(:,ik_ibz) = cweight(:,ik_ibz)*tetra%ibz_multiplicity(ik_ibz)/tetra%nkbz/tetra%tetra_total(ik_ibz)
   end do
 case(2)
   cweight = cweight*tetra%vv
 end select

 call xmpi_sum(cweight, comm, ierr)

end subroutine htetra_weights_wvals_zinv
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

end module m_htetrahedron
!!***

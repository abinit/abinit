!!****m* ABINIT/m_ddb
!! NAME
!!  m_ddb
!!
!! FUNCTION
!!  This module contains the declaration of data types and methods
!!  used to handle the blocks of data in DDB files:
!!  blkval, nrm, qpt, flg, and associated dimensions
!!  Main entry point for client code that needs to read the DDB data.
!!
!! COPYRIGHT
!! Copyright (C) 2011-2022 ABINIT group (MJV, XG, MT, MM, MVeithen, MG, PB, JCC, SP, GA)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_ddb

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
 use m_ddb_hdr
 use m_dtset
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use m_io_tools,       only : iomode_from_fname
 use defs_datatypes,   only : pseudopotential_type
 use m_fstrings,       only : sjoin, itoa, ktoa, endswith
 use m_numeric_tools,  only : mkherm
 use m_symtk,          only : mati3inv, matr3inv, littlegroup_q, symatm
 use m_io_tools,       only : get_unit
 use m_copy,           only : alloc_copy
 use m_geometry,       only : phdispl_cart2red, mkrdim, xred2xcart, metric
 use m_crystal,        only : crystal_t, crystal_init
 use m_dynmat,         only : cart29, d2sym3, cart39, d3sym, chneu9, asria_calc, asria_corr, asrprs, dfpt_phfrq, sytens
 use m_pawtab,         only : pawtab_type, pawtab_nullify, pawtab_free
 use m_psps,           only : psps_copy, psps_free

 implicit none

 private

 public :: rdddb9           ! This routine reads the derivative database entirely,
 public :: nlopt            ! Output of all quantities related to third-order derivatives of the energy.
 public :: chkin9
 public :: carttransf       ! Transform a second-derivative matrix (EIG2D) from reduced
                            ! coordinates to cartesian coordinates.
 public :: lwcart           ! Transform a 3rd order derivative tensor (long-wave) from reduced (actually
                            ! mixed since strain derivatives are already in cartesian) to cartesian
                            ! coordinates
 public :: ddb_lw_copy      ! Copy the ddb object after reading the long wave 3rd order derivatives
                            ! into a new ddb_lw and resizes ddb as for 2nd order derivatives

 public :: symdm9

 real(dp),public,parameter :: DDB_QTOL=2.0d-8
 ! Tolerance for the identification of two wavevectors
!!***

!----------------------------------------------------------------------

!!****t* m_ddb/ddb_type
!! NAME
!! ddb_type
!!
!! FUNCTION
!!  Provides methods to extract and post-process the results in the derivative database (DDB)
!!
!! SOURCE

 type,public :: ddb_type

  logical :: has_ncid_open
  ! Is currently reading a netcdf file

  integer :: iblock_d2eig_nc
  ! Is currently reading a netcdf file

  integer :: msize
  ! Maximum size of dynamical matrices and other perturbations (ddk, dde...)

  integer :: mpert
  ! Maximum number of perturbations

  integer :: nblok
  ! Number of 2dte blocks in present object

  integer :: natom
  ! Number of atoms in the unit cell.

  integer :: ntypat
  ! Number of type of atoms.

  integer :: occopt
  ! Occupation option.

  integer :: prtvol
  ! Verbosity level.

  integer :: nband
  ! Number of bands for eigenvalues derivatives
  ! This corresponds to d2eig arrary shape,
  ! but the actual number of band is nband / nsppol

  integer :: nkpt
  ! Number of k-points for eigenvalues derivatives

  ! GA: FIXME
  integer :: nsppol
  ! Number of spin components for eigenvalues derivatives
  ! This index is absorbed into nband, to limit array ranks to 7.

  integer :: current_iblok
  ! Number of k-points for eigenvalues derivatives

  ! These values are used to call the anaddb routines that don't use rprimd, gprimd.
  real(dp) :: rprim(3,3)
  real(dp) :: gprim(3,3)
  real(dp) :: acell(3)

  ! Many of these variables should become private so that one can refactor the ddb_t implementation
  integer,allocatable :: flg(:,:)
  ! flg(msize,nblok)
  ! Flag to indicate presence of a given block

  integer,allocatable :: typ(:)
  ! typ(nblok)
  ! Type of each block - nth-order derivatives of energy or eigenvalues.
  !      (0 => total energy)
  !      (1=> non-stationary block),
  !      (2=> stationary block),
  !      (3=> third order derivative).
  !      (4 => first-order derivatives of total energy)
  !      (5 => 2nd-order derivatives of eigenvalues)
  !      (33 => long wave third order derivatives of total energy)

  real(dp),allocatable :: amu(:)
  ! amu(ntypat)
  ! Mass of the atoms (atomic mass unit)

  real(dp),allocatable :: qpt(:,:)
  ! qpt(9,nblok)
  ! q-point vector in reciprocal space (reduced lattice coordinates) for each block
  ! Three possible phonon wavevectors can be specified for 3rd order derivatives,
  ! but only one should be used in case of second derivative of total energy,
  ! because we know that the second is the opposite of this value.

  real(dp),allocatable :: nrm(:,:)
  ! nrm(3,nblok)
  ! Normalization factors of the wavevectors for each block - can be 0 to indicate a direction of approach to gamma

  real(dp),allocatable :: val(:,:,:)
  ! val(2,msize,nblok)
  ! Values of the second energy derivatives in each block

  real(dp),allocatable :: kpt(:,:)
  ! kpt(3,nkpt)
  ! k-point vector in reciprocal space for eigenvalues derivatives

  real(dp),allocatable :: eig2dval(:,:,:,:)
  ! eig2dval(2,msize,nband,nkpt)
  ! Values of the second derivatives of eigenvalues
  ! Only a single block (a single q-point) is held in memory.
  ! Note that isppol index is wrapped into nband index.

  contains

    procedure :: init => ddb_init
     ! Construct the object from the dtset.

    procedure :: free => ddb_free
     ! Free dynamic memory.

    procedure :: malloc => ddb_malloc
     ! Allocate dynamic memory

    procedure :: malloc_d2eig => ddb_malloc_d2eig
     ! Allocate dynamic memory

    procedure :: copy => ddb_copy
     ! Copy the object.

    procedure :: set_qpt => ddb_set_qpt
     ! Set the wavevector

    procedure :: set_d1matr => ddb_set_d1matr
     ! Set values for the first-order derivative matrix in tensor shape

    procedure :: get_d1matr => ddb_get_d1matr
     ! Transform the first-order derivative matrix in tensor shape

    procedure :: set_d2matr => ddb_set_d2matr
     ! Set values for the second-order derivative matrix

    procedure :: get_d2matr => ddb_get_d2matr
     ! Transform the second-order derivative matrix in tensor shape

    procedure :: set_d3matr => ddb_set_d3matr
     ! Set values for the third-order derivative matrix

    procedure :: get_d3matr => ddb_get_d3matr
     ! Transform the third-order derivative matrix in tensor shape

    procedure :: get_d2eig => ddb_get_d2eig
     ! Transform the second-order derivative matrix of eigs in tensor shape

    procedure :: set_d2eig => ddb_set_d2eig
     ! Set values for the second-order derivative matrix of eigs

    procedure :: set_d2eig_reshape => ddb_set_d2eig_reshape
     ! Set values for the second-order derivative matrix of eigs
     ! with band index before perturbation indices

    procedure :: set_gred => ddb_set_gred
     ! Set the gradient of total energy in reduced coordinates

    procedure :: set_pel => ddb_set_pel
     ! Set the electronic polarization

    procedure :: set_strten => ddb_set_strten
     ! Set the stress tensor

    procedure :: set_etotal => ddb_set_etotal
     ! Set the total energy

    procedure :: set_brav => ddb_set_brav
     ! Set the bravais lattice.

    procedure :: bcast => ddb_bcast
     ! Broadcast the object.

    procedure :: get_etotal => ddb_get_etotal
     ! Read the GS total energy.

    procedure :: get_gred => ddb_get_gred
     ! Get the gradient of total energy in reduced coordinates

    procedure :: get_pel => ddb_get_pel
     ! Get the electronic polarization

    procedure :: get_strten => ddb_get_strten
     ! Get the stress tensor

    procedure :: get_dielt_zeff => ddb_get_dielt_zeff
     ! Reads the Dielectric Tensor and the Effective Charges

    procedure :: get_dielt => ddb_get_dielt
     ! Reads the Dielectric Tensor

    procedure :: get_quadrupoles => ddb_get_quadrupoles
     ! Reads the Quadrupoles

    procedure :: get_dchidet => ddb_get_dchidet
     ! Reads the non-linear optical susceptibility tensor and the
     ! first-order change in the linear dielectric susceptibility

    procedure :: diagoq => ddb_diagoq
     ! Compute the phonon frequencies at the specified q-point by performing
     ! a direct diagonalizatin of the dynamical matrix.

    procedure :: get_asrq0 => ddb_get_asrq0
     ! Return object used to enforce the acoustic sum rule

    procedure :: symmetrize_and_transform => ddb_symmetrize_and_transform
     ! Symmetrize, transform cartesian coordinates, and add missing components

    procedure :: write_block_txt => ddb_write_block_txt
     ! Writes blocks of data in the DDB in text format.

    procedure :: write => ddb_write
     ! Write the DDB file in either txt or netcdf format.

    procedure :: write_txt => ddb_write_txt
     ! Write the body of the DDB text file.

    procedure :: write_nc => ddb_write_nc
     ! Write the netcdf file (DDB.nc).

    procedure :: read_block_txt => ddb_read_block_txt
     ! Read blocks of data in the DDB.

    procedure :: get_block => ddb_get_block
     ! Finds the block containing the derivatives of the total energy.

    procedure :: read_d2eig => ddb_read_d2eig
     ! Read the next DDB block containing 2nd order derivatives of eigenvalues.

    procedure :: read_d2eig_txt => ddb_read_d2eig_txt
     ! Read the next DDB block containing 2nd order derivatives of eigenvalues.

    procedure :: read_d2eig_nc => ddb_read_d2eig_nc
     ! Read the next DDB block containing 2nd order derivatives of eigenvalues.

    procedure :: write_d2eig => ddb_write_d2eig
     ! Read the current DDB block containing 2nd order derivatives of eigenvalues.

    procedure :: write_d2eig_txt => ddb_write_d2eig_txt
     ! Read the current DDB block containing 2nd order derivatives of eigenvalues.

    procedure :: write_d2eig_nc => ddb_write_d2eig_nc
     ! Write the current DDB block containing 2nd order derivatives of eigenvalues.

    procedure :: read_d0E_nc => ddb_read_d0E_nc
     ! Read the next DDB block containing 0th order derivatives of energy.

    procedure :: read_d1E_nc => ddb_read_d1E_nc
     ! Read the next DDB block containing 1st order derivatives of energy.

    procedure :: read_d2E_nc => ddb_read_d2E_nc
     ! Read the next DDB block containing 2nd order derivatives of energy.

    procedure :: read_d3E_nc => ddb_read_d3E_nc
     ! Read the next DDB block containing 3rd order derivatives of energy.

    procedure :: from_file => ddb_from_file
     ! Construct the object from the DDB file.

    procedure :: read_txt => ddb_read_txt
     ! Construct the object from the DDB file in text format.

    procedure :: read_nc => ddb_read_nc
     ! Construct the object from the DDB file in netcdf format.

    procedure :: can_merge_blocks => ddb_can_merge_blocks
     ! Tell if two blocks can be merged

    procedure :: merge_blocks => ddb_merge_blocks
     ! Merge a block of an other ddb to the current object.

 end type ddb_type

 public :: ddb_to_dtset             ! Transfer ddb_hdr to dtset datatype
 public :: merge_ddb                ! Read a list of ddb files and merge them into a single ddb object

!!***

!!****t* m_ddb/asr_t
!! NAME
!!  asr_t
!!
!! FUNCTION
!!  Object used to enforce the acoustic sum rule from the Dynamical matrix at Gamma.
!!  Wraps several approaches that can be activated via the `asr` option.
!!
!! SOURCE

 type,public :: asrq0_t

   integer :: iblok = 0
    ! Index of the Gamma block in the DDB.
    ! Set to 0 if no block was found. Client code can use this flag to understand
    ! if ASR can be enforced.

   integer :: asr
   ! Option for the application of the ASR (input variable).

   integer :: natom
    ! Number of atoms.

   real(dp),allocatable :: d2asr(:,:,:,:,:)
   ! d2asr,(2,3,natom,3,natom))
   ! In case the interatomic forces are not calculated, the
   ! ASR-correction (d2asr) has to be determined here from the Dynamical matrix at Gamma.

   ! singular, uinvers and vtinvers are allocated and used only if asr in [3,4]
   ! i.e. Rotational invariance for 1D and 0D systems. dims=3*natom*(3*natom-1)/2
   real(dp),allocatable :: singular(:)
   ! singular,(1:dims))

   real(dp),allocatable :: uinvers(:,:)
   ! uinvers,(1:dims,1:dims))

   real(dp),allocatable :: vtinvers(:,:)
   ! vtinvers,(1:dims,1:dims))

 contains

   procedure :: apply => asrq0_apply
    ! Impose the acoustic sum rule based on the q=0 block found in the DDB file.

   procedure :: free => asrq0_free
    ! Free memory

 end type asrq0_t
!!***

 ! TODO: We should use this constants instead of magic numbers!
 ! BTW: Using a different value for NOSTAT and STAT is a non-sense!
 ! They both are 2-th order derivatives of the total energy!

 ! Flags used to indentify the block type.
 !integer,private,parameter :: DDB_BLKTYPE_ETOT = 0         ! Total energy
 !integer,private,parameter :: DDB_BLKTYPE_2DE_NOSTAT = 1   ! Second order derivative of the energy (non-stationary expression)
 !integer,private,parameter :: DDB_BLKTYPE_2DE_STAT = 2     ! Second order derivative of the energy (stationary expression)
 !integer,private,parameter :: DDB_BLKTYPE_3DE = 3          ! Third order derivative of the energy
 !integer,private,parameter :: DDB_BLKTYPE_1DE = 4          ! First order derivative of the energy
 !integer,private,parameter :: DDB_BLKTYPE_2DEIG = 5        ! Second order derivative of the eigenvalues

CONTAINS  !===========================================================
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_init
!! NAME
!! ddb_init
!!
!! FUNCTION
!!  Initialize a new ddb object for the current calculation.
!!
!! INPUTS
!!   ddb=the new ddb object
!!   dtset=dtset object of the current calculation
!!   nblok=number of blocks
!!   mpert=maximum number of perturbations (atom displacements + electric field + ...)
!!   with_d0E=this ddb contains 0th order derivatives
!!   with_d1E=this ddb contains 1st order derivatives
!!   with_d2E=this ddb contains 2nd order derivatives
!!   with_d3E=this ddb contains 3rd order derivatives
!!   with_d2eig=this ddb contains 2nd order derivatives of eigenvalues
!!   mband='number of bands' dimension of the d2eig array.
!!         Should actually correspond to the maximum number of bands for one kpoint
!!         multiplied by the number of spin polarization (mband*nsppol).
!!   nkpt=number of kpoints
!!   kpt=reduced coordinates of kpoints
!!
!! SOURCE

subroutine ddb_init(ddb, dtset, nblok, mpert, &
                    mband, nkpt, kpt,&
                    with_d0E, with_d1E, with_d2E, with_d3E, with_d2eig)

!Arguments -------------------------------
 class(ddb_type),intent(inout) :: ddb
 type(dataset_type),intent(in) :: dtset
 integer,intent(in) :: nblok, mpert
 integer,intent(in),optional :: mband,nkpt
 real(dp),intent(in),optional :: kpt(:,:)
 logical,intent(in),optional :: with_d0E, with_d1E, with_d2E, with_d3E, with_d2eig

!Local variables -------------------------------
 integer :: msize_, ii, ikpt, iblok
 logical :: with_d0E_, with_d1E_, with_d2E_, with_d3E_, with_d2eig_

! ************************************************************************

 with_d0E_   = .false. ; if (present(with_d0E))   with_d0E_ = with_d0E
 with_d1E_   = .false. ; if (present(with_d1E))   with_d1E_ = with_d1E
 with_d2E_   = .false. ; if (present(with_d2E))   with_d2E_ = with_d2E
 with_d3E_   = .false. ; if (present(with_d3E))   with_d3E_ = with_d3E
 with_d2eig_ = .false. ; if (present(with_d2eig)) with_d2eig_ = with_d2eig

 msize_ = 0
 if (with_d0E_) msize_ = 1
 if (with_d1E_) msize_ = 3 * mpert
 if (with_d2E_ .or. with_d2eig_) msize_ = 3 * mpert * 3 * mpert
 if (with_d3E_) msize_ = 3 * mpert * 3 * mpert * 3 * mpert

 call ddb%malloc(msize_, nblok, dtset%natom, dtset%ntypat, mpert)

 ddb%occopt = dtset%occopt
 ddb%prtvol = dtset%prtvol

 ddb%rprim(:,:) = dtset%rprim_orig(1:3,1:3,1)
 ddb%acell(:) = dtset%acell_orig(1:3,1)

 call matr3inv(ddb%rprim, ddb%gprim)

 ddb%qpt(:,:) = zero
 ddb%nrm(:,:) = one
 do iblok = 1,ddb%nblok
   if (with_d0E_) then
     ddb%typ(iblok) = BLKTYP_d0E_xx
   else if (with_d1E_) then
     ddb%typ(iblok) = BLKTYP_d1E_xx
   else if (with_d2E_) then
     ddb%typ(iblok) = BLKTYP_d2E_ns
   else if (with_d3E_) then
     ddb%typ(iblok) = BLKTYP_d3E_xx
   else if (with_d2eig_) then
     ddb%typ(iblok) = BLKTYP_d2eig_re
   end if
 end do
 ddb%flg(:,:) = 0
 ddb%amu(:) = dtset%amu_orig(:,1)

 ddb%nsppol = dtset%nsppol

 if (present(mband)) then
   ddb%nband = mband
 else
   ddb%nband = dtset%mband * ddb%nsppol
 end if

 if (present(nkpt)) then
   ddb%nkpt = nkpt
 else
   ddb%nkpt = dtset%nkpt
 end if

 ! TODO: Allocate d2eig here instead of leaving it to the calling routine.
 if (with_d2eig_) then
    call ddb%malloc_d2eig(ddb%nband, ddb%nkpt)
 end if

 if (present(kpt)) then
   do ikpt=1,ddb%nkpt
     do ii = 1,3
       ddb%kpt(ii,ikpt) = kpt(ii,ikpt)
     end do
   end do
 end if

end subroutine ddb_init
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_free
!! NAME
!! ddb_free
!!
!! FUNCTION
!!  Clean and deallocate types for the ddb_type structure
!!
!! SOURCE

subroutine ddb_free(ddb)

!Arguments -------------------------------
 class(ddb_type),intent(inout) :: ddb

! ************************************************************************

 !integer
 ABI_SFREE(ddb%flg)
 ABI_SFREE(ddb%typ)

 ! real
 ABI_SFREE(ddb%amu)
 ABI_SFREE(ddb%qpt)
 ABI_SFREE(ddb%nrm)
 ABI_SFREE(ddb%kpt)
 ABI_SFREE(ddb%val)
 ABI_SFREE(ddb%eig2dval)

end subroutine ddb_free
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_copy
!! NAME
!! ddb_copy
!!
!! FUNCTION
!!  Create object and copy all types for the ddb_type structure
!!
!! SOURCE

subroutine ddb_copy(iddb, oddb)

!Arguments -------------------------------
!array
 class(ddb_type),intent(in) :: iddb
 type(ddb_type),intent(out) :: oddb

! ************************************************************************

 ! Copy dimensions and static variables.
 oddb%msize = iddb%msize
 oddb%mpert = iddb%mpert
 oddb%nblok = iddb%nblok
 oddb%natom = iddb%natom
 oddb%ntypat = iddb%ntypat
 oddb%occopt = iddb%occopt
 oddb%prtvol = iddb%prtvol

 oddb%rprim = iddb%rprim
 oddb%gprim = iddb%gprim
 oddb%acell = iddb%acell

 ! Allocate and copy the allocatable arrays.
 call alloc_copy(iddb%flg, oddb%flg)
 call alloc_copy(iddb%typ, oddb%typ)
 call alloc_copy(iddb%amu, oddb%amu)
 call alloc_copy(iddb%nrm, oddb%nrm)
 call alloc_copy(iddb%qpt, oddb%qpt)
 call alloc_copy(iddb%val, oddb%val)

end subroutine ddb_copy
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_malloc
!! NAME
!! ddb_malloc
!!
!! FUNCTION
!!  Allocate dynamic memory.
!!
!! INPUTS
!!   msize=maximum size of one block of the ddb
!!         (e.g. 3*mpert * 3*mpert)
!!   nblok=number of blocks in the ddb
!!   natom=number of atoms
!!   ntypat=number of atom types
!!   mpert=maximum number of perturbations
!!         (atom displacements + electric field + ...)
!!   nkpt=number of k-points. Optional, indicates the use of eig2d.
!!   nband='number of bands' dimension of the d2eig array.
!!         Should actually correspond to the maximum number of bands for one kpoint
!!         multiplied by the number of spin polarization (mband*nsppol).
!!
!! OUTPUT
!!
!! SOURCE

subroutine ddb_malloc(ddb, msize, nblok, natom, ntypat, mpert, nkpt, nband)

!Arguments -------------------------------
!array
 class(ddb_type),intent(inout) :: ddb
 integer,intent(in) :: msize,nblok,natom,ntypat,mpert
 integer,intent(in),optional :: nkpt,nband

! ************************************************************************

 ddb%msize = msize
 ddb%nblok = nblok
 ddb%natom = natom
 !ddb%mpert = natom + MPERT_MAX
 ddb%mpert = mpert
 ddb%ntypat = ntypat

 ! integer
 ABI_CALLOC(ddb%flg, (msize, nblok))
 ABI_CALLOC(ddb%typ, (nblok))

 ! real
 ABI_MALLOC(ddb%amu, (ntypat))
 ABI_MALLOC(ddb%nrm, (3, nblok))
 ABI_MALLOC(ddb%qpt, (9, nblok))
 ABI_MALLOC(ddb%val, (2, msize, nblok))
 ddb%val = huge(one)

 ! FIXME: should really add nsppol argument (see thmeig).
 if (present(nkpt) .and. present(nband)) then
   call ddb%malloc_d2eig(nband, nkpt)
 end if

end subroutine ddb_malloc
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_malloc_d2eig
!! NAME
!! ddb_malloc_d2eig
!!
!! FUNCTION
!!  Allocate dynamic memory for second derivatives of eigenvalues.
!!
!! INPUTS
!!   mband='number of bands' dimension of the d2eig array.
!!         Should actually correspond to the maximum number of bands for one kpoint
!!         multiplied by the number of spin polarization (mband*nsppol).
!!   nkpt=number of kpoints
!!
!! OUTPUT
!!
!! SOURCE

subroutine ddb_malloc_d2eig(ddb, mband, nkpt)

!Arguments -------------------------------
!array
 class(ddb_type),intent(inout) :: ddb
 integer,intent(in) :: mband, nkpt

! ************************************************************************

  ddb%nband = mband / ddb%nsppol
  ddb%nkpt = nkpt
  ABI_MALLOC(ddb%kpt, (3, nkpt))
  ABI_MALLOC(ddb%eig2dval, (2, ddb%msize, mband, nkpt))

end subroutine ddb_malloc_d2eig
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_set_qpt
!! NAME
!! ddb_set_qpt
!!
!! FUNCTION
!!  Set the q-point wavevector for a certain block.
!!  In case of 3rd order derivatives, three q-points need to be specified
!!  with the constrain q1 + q2 + q3 = 0.
!!
!! INPUTS
!!  iblok=index of the block being set.
!!  qpt=reduced coordinates of first qpoint
!!  qpt2=reduced coordinates of second qpoint
!!  qpt3=reduced coordinates of third qpoint
!!
!! OUTPUT
!!
!! SOURCE

subroutine ddb_set_qpt(ddb, iblok, qpt, qpt2, qpt3)

!Arguments -------------------------------
!array
 class(ddb_type),intent(inout) :: ddb
 real(dp), intent(in) :: qpt(3)
 real(dp), intent(in),optional :: qpt2(3)
 real(dp), intent(in),optional :: qpt3(3)
!scalars
 integer,intent(in) :: iblok

! ************************************************************************

 ddb%qpt(1:3,iblok) = qpt(1:3)
 ddb%nrm(1,iblok) = one

 if (present(qpt2)) then
   ddb%qpt(4:6,iblok) = qpt2(1:3)
   ddb%nrm(2,iblok) = one
 end if

 if (present(qpt3)) then
   ddb%qpt(7:9,iblok) = qpt3(1:3)
   ddb%nrm(3,iblok) = one
 end if

end subroutine ddb_set_qpt
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_set_d2matr
!! NAME
!! ddb_set_d2matr
!!
!! FUNCTION
!!  Set values for the second-order derivative matrix.
!!
!! INPUTS
!!  iblok=index of the block being set.
!!  d2matr=the second-order derivative matrix.
!!  flg=flag to indicate presence of a given element.
!!
!! OUTPUT
!!
!! SOURCE

subroutine ddb_set_d2matr(ddb, iblok, d2matr, flg)

!Arguments -------------------------------
 class(ddb_type),intent(inout) :: ddb
!scalars
 integer,intent(in) :: iblok
!array
 real(dp), intent(in) :: d2matr(2,3,ddb%mpert,3,ddb%mpert)
 integer, intent(in) :: flg(3,ddb%mpert,3,ddb%mpert)

!Local variables -------------------------
!scalars
 integer :: idir1,idir2,ii,ipert1,ipert2

! ************************************************************************

 ii=0
 do ipert2=1,ddb%mpert
   do idir2=1,3
     do ipert1=1,ddb%mpert
       do idir1=1,3
         ii=ii+1
         ddb%flg(ii,iblok) = flg(idir1,ipert1,idir2,ipert2)
         ddb%val(1,ii,iblok) = d2matr(1,idir1,ipert1,idir2,ipert2)
         ddb%val(2,ii,iblok) = d2matr(2,idir1,ipert1,idir2,ipert2)
       end do
     end do
   end do
 end do

end subroutine ddb_set_d2matr
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_get_d2matr
!! NAME
!! ddb_get_d2matr
!!
!! FUNCTION
!!  Transform the second-order derivative matrix
!!  from flat indices to real tensor d2matr(cplex,ncart,natom,ncart,natom)
!!
!! INPUTS
!!  iblok=index of the block to get.
!!
!! OUTPUT
!!  d2matr=the second-order derivative matrix.
!!  flg=flag to indicate presence of a given element.
!!
!! SOURCE

subroutine ddb_get_d2matr(ddb, iblok, d2matr, flg)

!Arguments -------------------------------
!array
 class(ddb_type),intent(inout) :: ddb
 integer,intent(in) :: iblok
 real(dp), allocatable, intent(out) :: d2matr(:,:,:,:,:)
 integer, allocatable, intent(out) :: flg(:,:,:,:)
!scalars

!Local variables -------------------------
!scalars
 integer :: ii,idir1,idir2,ipert1,ipert2

! ************************************************************************

 ABI_MALLOC(d2matr, (2,3,ddb%mpert,3,ddb%mpert))
 ABI_MALLOC(flg, (3,ddb%mpert,3,ddb%mpert))

 d2matr = zero

 ii=0
 do ipert2=1,ddb%mpert
   do idir2=1,3
     do ipert1=1,ddb%mpert
       do idir1=1,3
         ii=ii+1
         flg(idir1,ipert1,idir2,ipert2) = ddb%flg(ii,iblok)
         if (ddb%flg(ii,iblok) > 0) then
           d2matr(1,idir1,ipert1,idir2,ipert2) = ddb%val(1,ii,iblok)
           d2matr(2,idir1,ipert1,idir2,ipert2) = ddb%val(2,ii,iblok)
         end if
       end do
     end do
   end do
 end do

end subroutine ddb_get_d2matr
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_set_gred
!! NAME
!! ddb_set_gred
!!
!! FUNCTION
!!  Set the forces in reduced coordinates (Hartree).
!!
!! INPUTS
!!  gred=the gradient of the total energy with respect
!!       to change of reduced coordinates
!!  iblok=index of the block being set.
!!
!! OUTPUT
!!
!! SOURCE

subroutine ddb_set_gred(ddb, gred, iblok)

!Arguments -------------------------------
!array
 class(ddb_type),intent(inout) :: ddb
 real(dp), intent(in) :: gred(3,ddb%natom)
!scalars
 integer,intent(in) :: iblok

!Local variables -------------------------
!scalars
 integer :: idir, iatom, indx

! ************************************************************************

 ddb%typ(iblok) = BLKTYP_d1E_xx
 indx = 0
 do iatom = 1, ddb%natom
   do idir = 1, 3
     indx = indx + 1
     ddb%flg(indx,iblok) = 1
     ddb%val(1,indx,iblok) = gred(idir,iatom)
     ddb%val(2,indx,iblok) = zero
   end do
 end do

end subroutine ddb_set_gred
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_set_pel
!! NAME
!! ddb_set_gred
!!
!! FUNCTION
!!  Set the electronic polarization.
!!
!! INPUTS
!!  pel=ucvol times the electronic polarization in reduced coordinates.
!!  flg=flag to indicate presence of a given element.
!!  iblok=index of the block being set.
!!
!! OUTPUT
!!
!! SOURCE

subroutine ddb_set_pel(ddb, pel, flg, iblok)

!Arguments -------------------------------
!array
 class(ddb_type),intent(inout) :: ddb
 real(dp), intent(in) :: pel(3)
 integer,intent(in) :: flg(3)
!scalars
 integer,intent(in) :: iblok

!Local variables -------------------------
!scalars
 integer :: idir, indx

! ************************************************************************

 ddb%typ(iblok) = BLKTYP_d1E_xx
 indx = 3*ddb%natom + 3
 do idir = 1, 3
   indx = indx + 1
   ddb%flg(indx,iblok) = flg(idir)
   ddb%val(1,indx,iblok) = pel(idir)
   ddb%val(2,indx,iblok) = zero
 end do

end subroutine ddb_set_pel
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_set_strten
!! NAME
!! ddb_set_strten
!!
!! FUNCTION
!!  Set the stress tensor.
!!
!! INPUTS
!!  strten=the stress tensor in cartesian coordinates.
!!  iblok=index of the block we are setting.
!!
!! OUTPUT
!!
!! SOURCE

subroutine ddb_set_strten(ddb, strten, iblok)

!Arguments -------------------------------
!array
 class(ddb_type),intent(inout) :: ddb
 real(dp), intent(in) :: strten(6)
!scalars
 integer,intent(in) :: iblok

!Local variables -------------------------
!scalars
 integer :: indx

! ************************************************************************

 ddb%typ(iblok) = BLKTYP_d1E_xx
 indx = 3*ddb%natom + 6

 ddb%flg(indx+1:indx+6,1) = 1
 ddb%val(1,indx+1:indx+6,1) = strten(1:6)
 ddb%val(2,indx+1:indx+6,1) = zero

end subroutine ddb_set_strten
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_get_d1matr
!! NAME
!! ddb_get_d1matr
!!
!! FUNCTION
!!  Transform the first-order derivative matrix
!!  from flat indices to real tensor d1matr(cplex,ncart,natom)
!!
!! INPUTS
!!  iblok=index of the block to get.
!!
!! OUTPUT
!!  d1matr=the first-order derivative matrix.
!!  flg=flag to indicate presence of a given element.
!!
!! SOURCE

subroutine ddb_get_d1matr(ddb, iblok, d1matr, flg)

!Arguments -------------------------------
!array
 class(ddb_type),intent(inout) :: ddb
 integer,intent(in) :: iblok
 real(dp), allocatable, intent(out) :: d1matr(:,:,:)
 integer, allocatable, intent(out) :: flg(:,:)
!scalars

!Local variables -------------------------
!scalars
 integer :: ii,idir1,ipert1

! ************************************************************************

 ABI_MALLOC(d1matr, (2,3,ddb%mpert))
 ABI_MALLOC(flg, (3,ddb%mpert))

 d1matr = zero

 ii=0
 do ipert1=1,ddb%mpert
   do idir1=1,3
     ii=ii+1
     flg(idir1,ipert1) = ddb%flg(ii,iblok)
     if (ddb%flg(ii,iblok) > 0) then
       d1matr(1,idir1,ipert1) = ddb%val(1,ii,iblok)
       d1matr(2,idir1,ipert1) = ddb%val(2,ii,iblok)
     end if
   end do
 end do

end subroutine ddb_get_d1matr
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_set_d1matr
!! NAME
!! ddb_set_d1matr
!!
!! FUNCTION
!!  Set values for the first-order derivative matrix.
!!
!! INPUTS
!!  iblok=index of the block being set.
!!  d1matr=the first-order derivative matrix.
!!  flg=flag to indicate presence of a given element.
!!
!! SOURCE

subroutine ddb_set_d1matr(ddb, iblok, d1matr, flg)

!Arguments -------------------------------
!array
 class(ddb_type),intent(inout) :: ddb
 real(dp), intent(in) :: d1matr(2,3,ddb%mpert)
 integer, intent(in) :: flg(3,ddb%mpert)
!scalars
 integer,intent(in) :: iblok

!Local variables -------------------------
!scalars
 integer :: ii,ipert1,idir1

! ************************************************************************

 ii=0
 do ipert1=1,ddb%mpert
   do idir1=1,3
     ii=ii+1
     ddb%val(1,ii,iblok) = d1matr(1,idir1,ipert1)
     ddb%val(2,ii,iblok) = d1matr(2,idir1,ipert1)
     ddb%flg(ii,iblok) = flg(idir1,ipert1)
   end do
 end do

end subroutine ddb_set_d1matr
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_set_etotal
!! NAME
!! ddb_set_etotal
!!
!! FUNCTION
!!  Set the total energy
!!
!! INPUTS
!!  etotal=the total energy.
!!  iblok=index of the block we are setting.
!!
!! OUTPUT
!!
!! SOURCE

subroutine ddb_set_etotal(ddb, etotal, iblok)

!Arguments -------------------------------
!array
 class(ddb_type),intent(inout) :: ddb
!scalars
 real(dp), intent(in) :: etotal
 integer,intent(in) :: iblok

! ************************************************************************

 ddb%typ(iblok) = BLKTYP_d0E_xx
 ddb%val(1,1,iblok) = etotal
 ddb%val(2,1,iblok) = zero
 ddb%flg(1,iblok) = 1

end subroutine ddb_set_etotal
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_set_brav
!! NAME
!! ddb_set_brav
!!
!! FUNCTION
!!  Modify the current values of rprim according to bravais lattice.
!!  Perform some checks on the primitive vectors
!!  before rescaling them such that rprim(1,2)=0.5
!!
!! INPUTS
!!  brav
!!   1 -> No rescaling.
!!   other -> Check and rescale.
!!
!!  Note that the meaning of brav is
!!    1 or -1 -> simple lattice
!!    2 -> face-centered cubic
!!    3 -> body-centered lattice
!!    4 -> hexagonal lattice (D6h)
!!
!! OUTPUT
!!
!! SOURCE

subroutine ddb_set_brav(ddb, brav)

!Arguments -------------------------------
!array
 class(ddb_type),intent(inout) :: ddb
!scalars
 integer,intent(in) :: brav

!Local variables-------------------------------
!scalars
 real(dp) :: factor
 character(len=500) :: msg

! *************************************************************************

 ! Renormalize rprim to possibly satisfy the constraint abs(rprim(1,2))=half when abs(brav)/=1
 ! This section is needed to preserver the behaviour of the old implementation.
 if (abs(brav)/=1 .and. abs(abs(ddb%rprim(1,2))-half)>tol10) then
   if(abs(ddb%rprim(1,2))<tol6)then
     write(msg, '(a,i0,7a)' )&
      'The input DDB value of brav is ',brav,',',ch10,&
      'and the one of rprim(1,2) is zero.',ch10,&
      'These are incompatible',ch10,&
      'Action: check the value of brav and rprim(1,2) in your DDB.'
     ABI_ERROR(msg)
   end if
   factor = abs(ddb%rprim(1,2)) * two
   ddb%acell(:) = ddb%acell(:) * factor
   ddb%rprim(:,:) = ddb%rprim(:,:) / factor
   ddb%gprim(:,:) = ddb%gprim(:,:) * factor
 end if

end subroutine ddb_set_brav
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_bcast
!! NAME
!! ddb_bcast
!!
!! FUNCTION
!!  MPI broadcast all types for the ddb_type structure
!!
!! INPUTS
!!   comm=MPI communicator
!!
!! SIDE EFFECTS
!!   Ddb<type(ddb_type)>= Input if node is master, other nodes returns with a completely initialized instance.
!!
!! SOURCE

subroutine ddb_bcast(ddb, comm)

!Arguments -------------------------------
!array
 class(ddb_type),intent(inout) :: ddb
 integer, intent(in) :: comm

!Local variables-------------------------------
!scalars
 integer, parameter :: master=0
 integer :: ierr
! *************************************************************************

 if (xmpi_comm_size(comm) == 1) return

 DBG_ENTER("COLL")

 ! Transmit dimensions and static variables.
 call xmpi_bcast(ddb%nblok, master, comm, ierr)
 call xmpi_bcast(ddb%natom, master, comm, ierr)
 call xmpi_bcast(ddb%ntypat, master, comm, ierr)
 call xmpi_bcast(ddb%nsppol, master, comm, ierr)
 call xmpi_bcast(ddb%mpert, master, comm, ierr)
 call xmpi_bcast(ddb%msize, master, comm, ierr)

 call xmpi_bcast(ddb%occopt, master, comm, ierr)
 call xmpi_bcast(ddb%prtvol, master, comm, ierr)

 !real
 call xmpi_bcast(ddb%rprim, master, comm, ierr)
 call xmpi_bcast(ddb%gprim, master, comm, ierr)
 call xmpi_bcast(ddb%acell, master, comm, ierr)

 ! Allocate arrays on the other nodes.
 if (xmpi_comm_rank(comm) /= master) then
   call ddb%malloc(ddb%msize, ddb%nblok, ddb%natom, ddb%ntypat, ddb%mpert)
 end if

 call xmpi_bcast(ddb%flg, master, comm, ierr)
 call xmpi_bcast(ddb%typ, master, comm, ierr)
 call xmpi_bcast(ddb%amu, master, comm, ierr)
 call xmpi_bcast(ddb%nrm, master, comm, ierr)
 call xmpi_bcast(ddb%qpt, master, comm, ierr)
 call xmpi_bcast(ddb%val, master, comm, ierr)

 DBG_EXIT("COLL")

end subroutine ddb_bcast
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_get_block
!!
!! NAME
!! ddb_get_block
!!
!! FUNCTION
!! This routine finds the block that contains the
!! information on the derivatives of the total energy specified
!! by the parameters rfphon,rfelfd,rfstrs,rftyp and
!! the phonon wavevectors qphon (and their normalisation).
!! In case the DDB does not contain this information, the subroutine returns iblok=0
!!
!! INPUTS
!! ddb = ddb blok datastructure
!!   flg(msize,nblok)=flag for every matrix element:
!!      0 => the element is not in the data block.
!!      1 => the element is in the data blok.
!!   nrm(3,nblok)=normalization factors for the three allowed wavevectors
!!   qpt(3,nblok)=wavevector of the perturbation(s). The elements
!!   typ(nblok)=type of the block.
!!      (1=> non-stationary block),
!!      (2=> stationary block),
!!      (3=> third order derivative).
!! qphon(3,3)=wavevectors for the three possible phonons
!!  (note : only one should be used in case of second derivative of total energy,
!!  because we know that the second is the opposite of this value)
!! qphnrm(3) =normalisation factors for the three possible phonons
!! rfphon(4) = 1=> response to phonons
!!             2=> second derivative of total energy
!! rfelfd(4) = 1=> d/dk, 2=> electric field only, 3=> both (see comment on rfphon)
!! rfstrs(4) = 1=> uniaxial stresses, 2=> shear stresses, 3=> both (see comment on rfphon)
!! rftyp =
!!   0 => total energy
!!   1 => non-stationary formulation of the 2nd derivative
!!   2 => stationary formulation of the 2nd derivative
!!   3 => third derivative of total energy
!!   4 => first-order derivatives of total energy
!!  33 => long wave third order derivatives of total energy
!! [rfqvec(4)] = 1=> d/dq (optional)
!!
!! OUTPUT
!! iblok= number of the block that corresponds to the specifications. 0 if not found.
!!
!! SOURCE

subroutine ddb_get_block(ddb, iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, rftyp, rfqvec)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: rftyp
 integer,intent(out) :: iblok
 class(ddb_type),intent(in) :: ddb
!arrays
 integer,intent(in) :: rfelfd(4),rfphon(4),rfstrs(4)
 real(dp),intent(inout) :: qphnrm(3),qphon(3,3)
 integer,optional,intent(in) :: rfqvec(4)

!Local variables -------------------------
!scalars
 integer :: blkgam,ider,idir,idir1,idir2,idir3,ii,index,ipert,ipert1,ipert2
 integer :: ipert3,nder,ok,mpert,natom
 character(len=500) :: msg
!arrays
 integer :: gamma(3)
 integer,allocatable :: worki(:,:)
 real(dp) :: qpt(3)
 integer :: rfqvec_(4)

! *********************************************************************

 mpert = ddb%mpert
 natom = ddb%natom

 ! Get the number of derivative
 if (is_type_d2E(rftyp)) then
   nder=2
 else if (is_type_d3E(rftyp)) then
   nder=3
 else if (is_type_d0E(rftyp)) then
   nder=0
 else if (is_type_d1E(rftyp)) then
   nder=1
 else
   write(msg, '(a,i0,a)')' rftyp is equal to ',rftyp,'. The only allowed values are 0, 1, 2, 3 or 4.'
   ABI_BUG(msg)
 end if

 rfqvec_(:)=0; if(present(rfqvec))rfqvec_(:)=rfqvec(:)

 ! In case of a second-derivative, a second phonon wavevector is provided.
 if(nder==2)then
   do ii=1,3
     qphon(ii,2)=-qphon(ii,1)
   end do
   qphnrm(2)=qphnrm(1)
 end if

 ! In case of a third derivative, the sum of wavevectors to gamma is checked
 if (nder == 3) then
   qpt(:) = qphon(:,1)/qphnrm(1) + qphon(:,2)/qphnrm(2) + qphon(:,3)/qphnrm(3)
   call gamma9(gamma(nder),qpt,qphnrm(1),DDB_QTOL)
   if (gamma(nder) == 0) then
     write(msg,'(a,a,a)')&
      'the sum of the wavevectors of the third-order energy is ',ch10,&
      'not equal to zero'
     ABI_ERROR(msg)
   end if
 end if

 ! Check the validity of the requirement
 do ider=1,nder
   ! Identifies if qphon is at gamma
   call gamma9(gamma(ider),qphon(1:3,ider),qphnrm(ider),DDB_QTOL)

   if(gamma(ider)==0)then
     if(rfstrs(ider)/=0.or.rfelfd(ider)/=0.or.rfqvec_(ider)/=0)then
       write(msg, '(a,a)' )&
        'Not yet able to handle stresses or electric fields',ch10,&
        'with non-zero wavevector.'
       ABI_BUG(msg)
     end if
   end if
 end do

 ! Initialise the perturbation table
 ABI_MALLOC(worki,(mpert,4))
 worki(:,1:nder)=0

 ! Build the perturbation table
 do ider=1,nder
   ! First the phonons
   if(rfphon(ider)==1)then
     do ipert=1,natom
       worki(ipert,ider)=1
     end do
   end if
   ! Then the d/dk
   if (rfelfd(ider)==1.or.rfelfd(ider)==3) worki(natom+1,ider)=1
   ! Then the electric field
   if (rfelfd(ider)==2.or.rfelfd(ider)==3) worki(natom+2,ider)=1
   ! Then the ddq
   if (rfqvec_(ider)==1) worki(natom+8,ider)=1
   ! Then the uniaxial stress
   if (rfstrs(ider)==1.or.rfstrs(ider)==3) worki(natom+3,ider)=1
   ! At last, the shear stress
   if(rfstrs(ider)==2.or.rfstrs(ider)==3) worki(natom+4,ider)=1
 end do

 ! Examine every blok:
 do iblok=1,ddb%nblok

   ! If this variable is still 1 at the end of the examination, the blok is the good one...
   ok=1

   ! Check the type
   if(rftyp/=ddb%typ(iblok)) ok=0

   ! Check the wavevector
   if( ok==1 )then

     if (nder == 2) then
       call gamma9(blkgam,ddb%qpt(1:3,iblok),ddb%nrm(1,iblok),DDB_QTOL)
       if(blkgam/=gamma(1))then
         ok=0
       else if(blkgam==0)then
         do idir=1,3
           if( abs( ddb%qpt(idir,iblok)/ddb%nrm(1,iblok) - qphon(idir,1)/qphnrm(1) )>DDB_QTOL ) ok=0
         end do
       end if

     else if (nder == 3) then
       do ider = 1, nder
         do idir=1,3
           if( abs( ddb%qpt(idir+3*(ider-1),iblok)/ddb%nrm(ider,iblok) - qphon(idir,ider)/qphnrm(ider) )>DDB_QTOL )then
             ok=0
           end if ! qphon
         end do ! idir
       end do ! nder
     end if  ! nder

   end if ! ok

   ! Check if there is enough information in this blok
   if( ok==1 )then

     if (nder == 0) then
       if (ddb%flg(1,iblok) /= 1) then
         ok = 0
         if (ddb%prtvol > 1) then
           write(msg,'(a,i0,3a)' )&
            'The block ',iblok,' does not match the requirement',ch10,&
            'because it lacks the total energy'
           ABI_COMMENT(msg)
         end if
       end if
     end if

     do ipert1=1,mpert

       if ((nder == 4).and.(worki(ipert1,4) == 1).and.(ok == 1)) then
         do idir1 = 1, 3
           index = 3*(ipert1 - 1) + idir1
           if (ddb%flg(index,iblok) /= 1) ok = 0
         end do
       end if

       if (worki(ipert1,1)==1 .and. ok==1 )then
         do ipert2=1,mpert
           if (worki(ipert2,2)==1 .and. ok==1 )then
             do idir1=1,3
               do idir2=1,3

                 if (nder == 2) then
                   index=idir1+ 3*((ipert1-1)+mpert*((idir2-1)+3*(ipert2-1)))
                   if (ddb%flg(index,iblok)/=1) ok=0

                 else if (nder == 3) then
                   do ipert3 = 1, mpert
                     if (worki(ipert3,3) == 1 .and. ok == 1) then
                       do idir3 = 1, 3
                         index = idir1 + &
                           3*((ipert1 - 1) + mpert*((idir2 - 1) + &
                           3*((ipert2 -1 ) + mpert*((idir3 - 1) + 3*(ipert3 - 1)))))
                         if (ddb%flg(index,iblok) /= 1) ok = 0
                       end do  ! idir3
                     end if ! worki(ipert3,3)
                   end do ! i3pert
                 end if

               end do
             end do
           end if
         end do
       end if
     end do
   end if

   ! Now that everything has been checked, eventually end the search
   if(ok==1)exit
 end do

 if(ok==0)then
   iblok=0

   if (ddb%prtvol > 1) then
     write(msg, '(3a)' )&
      ' gtblk9 : ',ch10,&
      '  Unable to find block corresponding to the following specifications :'
     call wrtout(std_out,msg)
     write(msg, '(a,i3)' )' Type (rfmeth) =',rftyp
     call wrtout(std_out,msg)
     write(msg, '(a)' ) ' ider qphon(3)         qphnrm   rfphon rfelfd rfstrs rfqvec'
     call wrtout(std_out,msg)
     do ider=1,nder
       write(msg, '(i4,4f6.2,4i7)' )&
       ider,(qphon(ii,ider),ii=1,3),qphnrm(ider),rfphon(ider),rfelfd(ider),rfstrs(ider),rfqvec_(ider)
       call wrtout(std_out,msg)
     end do
   end if
 end if

 if (ok==1 .and. ddb%prtvol > 1) then
   write(msg,'(a,i0,2a)')' gtblk9: found block number ',iblok,' agree with',' specifications '
   call wrtout(std_out,msg)
 end if

 ABI_FREE(worki)

end subroutine ddb_get_block
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/gamma9
!!
!! NAME
!! gamma9
!!
!! FUNCTION
!! This small routine checks if the wavevector qphon and the
!! corresponding normalisation factor represent a phonon at Gamma.
!!
!! INPUTS
!!  qphon(3)=wavevector
!!  qphnrm=normalisation factor
!!  qtol=tolerance
!!
!! OUTPUT
!! gamma= if 1, means that the wavevector is indeed at Gamma otherwise 0.
!!
!! SOURCE


subroutine gamma9(gamma,qphon,qphnrm,qtol)

!Arguments -------------------------------
!scalars
 integer,intent(out) :: gamma
 real(dp),intent(in) :: qphnrm,qtol
!arrays
 real(dp),intent(in) :: qphon(3)

! *********************************************************************

 if( (abs(qphon(1))<qtol .and. abs(qphon(2))<qtol .and. abs(qphon(3))<qtol) .or. abs(qphnrm)<qtol ) then
   gamma=1
 else
   gamma=0
 end if

end subroutine gamma9
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_read_block_txt
!!
!! NAME
!! ddb_read_block_txt
!!
!! FUNCTION
!! Read the next block of data from a DDB in text format.
!!
!! INPUTS
!!  iblok=the blok index to be assigned
!!  mpert=maximum number of ipert
!!  msize=maximum size of the arrays flags and values
!!  nunit=unit number for the data block file
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output
!! ddb = ddb block datastructure
!! ddb%typ=type of the block:
!!   0 => total energy
!!   1 => second-order energy derivatives, non-stationary block
!!   2 => second-order energy derivatives, stationary block
!!   3 => third-order energy derivatives
!!   4 => first-order energy derivatives: forces, stresses and polarization
!!   5 => second-order eigenvalue derivatives
!! ddb%flg(msize)=flag for every matrix element (0=> the element is
!!  not in the data block), (1=> the element is in the data blok)
!! ddb%qpt(9)=wavevector of the perturbation(s). The elements from
!!  1 to 3 are used if we are dealing with the 2nd derivative of
!!  total energy (only one wavevector), while all elements are
!!  used in case of a third order derivative of total energy (three wavevector could be present)
!! ddb%nrm(3)=normalization factors for the three allowed wavevectors.
!! ddb%val(2,msize)=real(dp), complex, value of the matrix elements that are present in the data block
!! [blkval2(2,msize,mband,nkpt)]= value of the matrix elements that are present in a block of EIGR2D/EIGI2D
!!
!! NOTES
!! only executed by one processor.
!!
!! SOURCE

subroutine ddb_read_block_txt(ddb,iblok,mband,mpert,msize,nkpt,nunit,&
                          blkval2,kpt) !optional

!Arguments -------------------------------
!scalars
 integer,intent(in) :: mband,mpert,msize,nkpt,nunit
 integer, intent(in) :: iblok
 !logical, intent(in), optional :: eig2d
 class(ddb_type),intent(inout) :: ddb
!arrays
 real(dp),intent(out),optional :: kpt(3,nkpt)
 real(dp),intent(out),optional :: blkval2(2,msize,mband,nkpt)

!Local variables -------------------------
!scalars
 integer :: band,iband,idir1,idir2,idir3,ii,ikpt,index,ipert1,ipert2,ipert3,nelmts
 logical :: eig2d_
 real(dp) :: ai,ar
 character(len=32) :: name
 character(len=500) :: msg

! *********************************************************************

 ! Zero every flag
 ddb%flg(1:msize, iblok)=0

 eig2d_ = .false.


 if(present(kpt).and.present(blkval2)) then
   ! GA: Weird that it is not allocated here
   blkval2(:,:,:,:)=zero
   kpt(:,:)=zero
   eig2d_ = .true.
 end if

 !if (present(eig2d)) then
 !   eig2d_ = eig2d
 !end if

 ! Read the block type and number of elements
 read(nunit,*)
 read(nunit, '(a32,12x,i8)' )name,nelmts

 ! TODO: Replace with STRING_d2E, etc.
 ! GA: Note that older versions used the expression '2rd' instead of '2nd'
 ! So this substitution needs to be checked for backward compatibility.
 ! Also, the strings for d2eig and d2eig_brd are undistinguishable
 ! I don't think the d2eig_brd was ever read by abinit or anaddb.
 if(name==' 2nd derivatives (non-stat.)  - ' .or. name==' 2rd derivatives (non-stat.)  - ')then
   ddb%typ(iblok)=BLKTYP_d2E_ns
 else if(name==' 2nd derivatives (stationary) - ' .or. name==' 2rd derivatives (stationary) - ')then
   ddb%typ(iblok)=BLKTYP_d2E_st
 else if(name==' 3rd derivatives              - ')then
   ddb%typ(iblok)=BLKTYP_d3E_xx
 else if(name==' Total energy                 - ')then
   ddb%typ(iblok)=BLKTYP_d0E_xx
 else if(name==' 1st derivatives              - ')then
   ddb%typ(iblok)=BLKTYP_d1E_xx
 else if(name==' 2nd eigenvalue derivatives   - ' .or. name==' 2rd eigenvalue derivatives   - ')then
   ddb%typ(iblok)=BLKTYP_d2eig_re
 else if(name==' 3rd derivatives (long wave)  - ')then
   ddb%typ(iblok)=BLKTYP_d3E_lw
 else
   write(msg,'(6a)')&
   'The following string appears in the DDB in place of',&
   ' the block type description :',ch10,trim(name),ch10,&
   'Action: check your DDB.'
   ABI_ERROR(msg)
 end if

 ! Read the 2nd derivative block
 if (is_type_d2E(ddb%typ(iblok))) then

   ! First check if there is enough space to read it
   if(msize<(3*mpert*3*mpert))then
     write(msg,'(3a)')&
     'There is not enough space to read a second-derivative block.',ch10,&
     'Action: increase msize and recompile.'
     ABI_ERROR(msg)
   end if

   ! Read the phonon wavevector
   read(nunit, '(4x,3es16.8,f6.1)' )(ddb%qpt(ii,iblok),ii=1,3),ddb%nrm(1,iblok)

   ! Read every element
   do ii=1,nelmts
     read(nunit,*)idir1,ipert1,idir2,ipert2,ar,ai
     index=idir1+3*((ipert1-1)+mpert*((idir2-1)+3*(ipert2-1)))
     ddb%flg(index,iblok)=1
     ddb%val(1,index,iblok)=ar
     ddb%val(2,index,iblok)=ai
   end do

 else if (is_type_d3E(ddb%typ(iblok))) then
   ! Read the 3rd derivative block

   ! First check if there is enough space to read it
   if(msize<(3*mpert*3*mpert*3*mpert))then
     write(msg, '(a,a,a,i10,a,i10,a,a,a)' )&
     'There is not enough space to read a third-derivative block.',ch10,&
     'The size provided is only ',msize,' although ',3*mpert*3*mpert*3*mpert,' is needed.',ch10,&
     'Action: increase msize and recompile.'
     ABI_ERROR(msg)
   end if

   ! Read the perturbation wavevectors
   read(nunit,'(4x,3es16.8,f6.1)')(ddb%qpt(ii,iblok),ii=1,3),ddb%nrm(1,iblok)
   read(nunit,'(4x,3es16.8,f6.1)')(ddb%qpt(ii,iblok),ii=4,6),ddb%nrm(2,iblok)
   read(nunit,'(4x,3es16.8,f6.1)')(ddb%qpt(ii,iblok),ii=7,9),ddb%nrm(3,iblok)

   ! Read every element
   do ii=1,nelmts
     read(nunit,'(6i4,2d22.14)')idir1,ipert1,idir2,ipert2,idir3,ipert3,ar,ai
     index=idir1+                     &
       3*((ipert1-1)+mpert*((idir2-1)+ &
       3*((ipert2-1)+mpert*((idir3-1)+3*(ipert3-1)))))
     ddb%flg(index,iblok)=1
     ddb%val(1,index,iblok)=ar
     ddb%val(2,index,iblok)=ai
   end do


 else if (is_type_d0E(ddb%typ(iblok))) then
   ! Read the total energy
   ! First check if there is enough space to read it
   if(msize<1)then
     write(msg, '(3a,i0,3a)' )&
      'There is not enough space to read a total energy block.',ch10,&
      'The size provided is only ',msize,' although 1 is needed.',ch10,&
      'Action: increase msize and recompile.'
     ABI_ERROR(msg)
   end if

   ! Read the total energy
   read(nunit,'(2d22.14)')ar,ai
   ddb%flg(1,iblok)=1
   ddb%val(1,1,iblok)=ar
   ddb%val(2,1,iblok)=ai


 else if (is_type_d1E(ddb%typ(iblok))) then
   !  Read the 1st derivative block
   !  First check if there is enough space to read it
   if (msize < (3*mpert)) then
     write(msg, '(3a,i0,a,i0,3a)' )&
     'There is not enough space to read a first-derivative block.',ch10,&
     'The size provided is only ',msize,' although ',3*mpert,' is needed.',ch10,&
     'Action: increase msize and recompile.'
     ABI_ERROR(msg)
   end if

   ! Read every element
   do ii=1,nelmts
     read(nunit,'(2i4,2d22.14)')idir1,ipert1,ar,ai
     index=idir1 + 3*(ipert1 - 1)
     ddb%flg(index,iblok)=1
     ddb%val(1,index,iblok)=ar
     ddb%val(2,index,iblok)=ai
   end do


 else if (is_type_d2eig(ddb%typ(iblok))) then

   ! Read the 2nd eigenvalue derivative block
   ! First check if there is enough space to read it
   if(msize<(3*mpert*3*mpert))then
     write(msg, '(3a,i0,a,i0,3a)' )&
     'There is not enough space to read a second-derivative block.',ch10,&
     'The size provided is only ',msize,' although ',3*mpert*3*mpert*mband*nkpt,' is needed.',ch10,&
     'Action: increase msize and recompile.'
     ABI_ERROR(msg)
   end if

   ! Read the phonon wavevector
   read(nunit, '(4x,3es16.8,f6.1)' )(ddb%qpt(ii,iblok),ii=1,3),ddb%nrm(1,iblok)

   ! Read the K point and band
   if (eig2d_) then
     do ikpt=1,nkpt
       read(nunit, '(9x,3es16.8)')(kpt(ii,ikpt),ii=1,3)
       do iband=1,mband
         read(nunit, '(6x,i3)') band
         ! Read every element
         do ii=1,nelmts
           read(nunit,*)idir1,ipert1,idir2,ipert2,ar,ai
           index=idir1+3*((ipert1-1)+mpert*((idir2-1)+3*(ipert2-1)))
           ddb%flg(index,iblok)=1
           blkval2(1,index,iband,ikpt)=ar
           blkval2(2,index,iband,ikpt)=ai
         end do !nelmts
       end do  !band
     end do !kpt
   end if
 end if

end subroutine ddb_read_block_txt
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_read_d2eig
!!
!! NAME
!! ddb_read_d2eig
!!
!! FUNCTION
!! Read the next DDB block containing second-order derivatives of eigenvalues
!! and store it in block number iblok.
!! The values of nband and nkpt must be set.
!!
!! INPUTS
!!  ddb_hdr=ddb header object with open file.
!!  iblok_store=the block index in the ddb object
!!  iblok_read=the block index in the ddb file
!!
!! OUTPUT
!!
!! NOTE
!! The ddb object must be allocated becore calling this routine.
!!
!! SOURCE


subroutine ddb_read_d2eig(ddb, ddb_hdr, iblok_store, iblok_read, comm)

!Arguments -------------------------------
!scalars
 class(ddb_type),intent(inout) :: ddb
 type(ddb_hdr_type),intent(in) :: ddb_hdr
 integer, intent(in) :: iblok_store
 integer, intent(in),optional :: iblok_read
 integer, intent(in),optional :: comm

!Local variables -------------------------
!scalars
 integer,parameter :: master=0
 integer :: comm_
 character(len=500) :: msg

! *********************************************************************

  if (present(comm)) then
    comm_ = comm
  else
    comm_ = xmpi_comm_self
  end if

  if (xmpi_comm_rank(comm_) == master) then

    if (ddb_hdr%has_open_file_nc) then

      ! Read the specified block and store it
      call ddb%read_d2eig_nc(ddb_hdr%ncid, iblok_store, iblok_read)

    else if (ddb_hdr%has_open_file_txt) then

      ! Read the next block and store it
      call ddb%read_d2eig_txt(ddb_hdr%unddb, iblok_store)

    else 
      write(msg, '(3a)' )&
      ! File has not beed open by ddb_hdr
      'Attempting to read from unopen file DDB.',ch10,&
      'Action: contact Abinit group.'
      ABI_ERROR(msg)
    end if

  end if

  ! GA: Should in principle broadcast the d2eig array,
  !     but I dont think it is required yet.

end subroutine ddb_read_d2eig
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_read_d2eig_txt
!!
!! NAME
!! ddb_read_d2eig_txt
!!
!! FUNCTION
!! Read the next DDB block containing second-order derivatives of eigenvalues
!! and store it in block number iblok.
!! The ddb object must have been allocated with nband and nkpt.
!!
!! INPUTS
!!  unddb=unit for the open file in text format
!!  iblok=the block index in the ddb object
!!
!! OUTPUT
!!
!! SOURCE


subroutine ddb_read_d2eig_txt(ddb, unddb, iblok)

!Arguments -------------------------------
!scalars
 class(ddb_type),intent(inout) :: ddb
 integer, intent(in) :: unddb
 integer, intent(in), optional :: iblok
!Local variables -------------------------
!scalars
 integer :: iblok_eig2d

! *********************************************************************

  iblok_eig2d = 1
  if (present(iblok)) iblok_eig2d = iblok

   ! GA: Here, nband should really be nband * nsppol.
   !     but this is the responsibility of the calling routine
   !     see thmeig and merge_ddb
   !     FIXME This is inconsistent with ddb_malloc_d2eig...
   !     I think I should change this with ddb%nband * ddb%nsppol
  call ddb%read_block_txt(iblok_eig2d,ddb%nband,ddb%mpert,ddb%msize,ddb%nkpt,unddb,&
                      ddb%eig2dval(:,:,:,:),ddb%kpt(:,:))

end subroutine ddb_read_d2eig_txt
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/rdddb9
!! NAME
!! rdddb9
!!
!! FUNCTION
!! This routine reads the derivative database entirely,
!! for use in ppddb9, and performs some checks and symmetrisation
!! At the end, the whole DDB is in central memory, contained in the array ddb%val(2,msize,ddb%nblok).
!!
!! The information on it is contained in the four arrays
!!   ddb%flg(msize,ddb%nblok) : blok flag for each element
!!   ddb%qpt(9,ddb%nblok)  : blok wavevector (unnormalized)
!!   ddb%nrm(3,ddb%nblok)  : blok wavevector normalization
!!   ddb%typ(ddb%nblok)    : blok type
!!
!! INPUTS
!! unddb = unit number for DDB io
!! dimekb=dimension of ekb (for the time being, only for norm- conserving psps)
!! iout=unit number for output of formatted data
!! lmnmax=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!       =if useylm=0, max number of (l,n)   comp. over all type of psps
!! mband=maximum number of bands
!! mpert =maximum number of ipert
!! msize=maximum size of data blocks
!! msym =maximum number of symmetry elements in space group
!! natom = number of atoms
!! ntypat=number of atom types
!! usepaw= 0 for non paw calculation; =1 for paw calculation
!! [raw] = 1 -> do not perform any symetrization or transformation to cartesian coordinates.
!!         0 (default) -> do perform these transformations.
!!
!! OUTPUT
!! acell(3)=length scales of cell (bohr)
!! amu(ntypat)=mass of the atoms (atomic mass unit)
!! ddb: ddb blok datatype
!!   contents: ddb%flg(msize,nblok)= flag of existence for each element of the DDB
!!             ddb%nrm(3,nblok)  : blok wavevector normalization
!!             ddb%qpt(9,nblok)  : blok wavevector (unnormalized)
!!             ddb%typ(nblok)    : blok type
!!             ddb%val(2,msize,nblok)= value of each complex element of the DDB
!!             ddb%nblok= number of bloks in the DDB
!! gmet(3,3)=reciprocal space metric tensor in bohr**-2
!! gprim(3,3)=dimensionless reciprocal space primitive translations
!! indsym(4,msym,natom)=indirect indexing array for symmetries
!! natom=number of atoms in cell
!! nsym=number of space group symmetries
!! rmet(3,3)=metric tensor in real space (bohr^2)
!! rprim(3,3)= primitive translation vectors
!! symrec(3,3,nsym)=3x3 matrices of the group symmetries (reciprocal space)
!! symrel(3,3,nsym)=3x3 matrices of the group symmetries (real space)
!! symafm(nsym)=Anti-ferromagnetic symmetries.
!! tnons(3,nsym)=fractional nonsymmorphic translations
!! typat(natom)=type integer for each atom in cell
!! ucvol=unit cell volume in bohr**3
!! xcart(3,natom)=atomic cartesian coordinates
!! xred(3,natom)=fractional dimensionless atomic coordinates
!! zion(ntypat)=charge on each type of atom (real number)
!! znucl(ntypat)=Nuclear charge for each type of pseudopotential
!!
!! SOURCE

subroutine rdddb9(ddb,ddb_hdr,unddb,&
                  acell,amu,gmet,gprim,indsym,&
                  mband,mpert,msize,msym,natom,nkpt,nsym,ntypat,&
                  rmet,rprim,symrec,symrel,symafm,tnons,typat,ucvol,&
                  xcart,xred,zion,znucl,raw)

!Arguments -------------------------------
! NOTE: these are used for dimensioning and then re-assigned in ioddb8.
!   This is almost definitely bad practice. In particular
!    it should be indsym(4,msym,natom),
!   and
!    the allocation allocate(kpt(3,nkpt)) is strange
!scalars
 integer,intent(in) :: unddb,mband,mpert,msize,msym
 integer,intent(inout) :: natom,nkpt,nsym,ntypat
 real(dp),intent(out) :: ucvol
 type(ddb_type),intent(inout) :: ddb
 type(ddb_hdr_type),intent(inout) :: ddb_hdr
 integer,optional,intent(in) :: raw
!arrays
 integer,intent(inout) :: indsym(4,msym,natom)
 integer,intent(out) :: symrec(3,3,msym),symrel(3,3,msym),symafm(msym)
 integer,intent(out) :: typat(natom)
 real(dp),intent(out) :: acell(3),amu(ntypat)
 real(dp),intent(out) :: gmet(3,3),gprim(3,3),rmet(3,3)
 real(dp),intent(out) :: rprim(3,3),tnons(3,msym),xcart(3,natom),xred(3,natom)
 real(dp),intent(out) :: zion(ntypat),znucl(ntypat)

!Local variables -------------------------
!mtyplo=maximum number of type, locally
!scalars
 integer,parameter :: msppol=2,mtyplo=6
 integer :: raw_
 integer :: iblok,isym
 real(dp),parameter :: tolsym8=tol8
!arrays
 real(dp) :: gprimd(3,3),rprimd(3,3)

! *********************************************************************

 DBG_ENTER("COLL")

 if (present(raw)) then
   raw_ = raw
 else
   raw_ = 0
 end if

 ! FIXME
 ! GA: Most of this stuff could be moved up to the calling routine

 nsym = ddb_hdr%nsym
 acell = ddb_hdr%acell
 rprim = ddb_hdr%rprim

 amu(:) = ddb_hdr%amu(1:ntypat)
 typat(:) = ddb_hdr%typat(1:natom)
 zion(:) = ddb_hdr%zion(1:ntypat)
 znucl(:) = ddb_hdr%znucl(1:ntypat)

 symafm(:) = ddb_hdr%symafm(:)
 symrel(:,:,:) = ddb_hdr%symrel(:,:,:)
 tnons(:,:) = ddb_hdr%tnons(:,:)

 xred(:,:) = ddb_hdr%xred(:,:)

 !call ddb_hdr%free()

 ! Compute different matrices in real and reciprocal space, also
 ! checks whether ucvol is positive.
 call mkrdim(acell,rprim,rprimd)

 ! call metric without printing to output
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 ! Obtain reciprocal space primitive transl g from inverse trans of r
 ! (Unlike in abinit, gprim is used throughout ifc; should be changed, later)
 call matr3inv(rprim,gprim)

 ! Generate atom positions in cartesian coordinates
 call xred2xcart(natom,rprimd,xcart,xred)

 ! Transposed inversion of the symmetry matrices, for use in the reciprocal space
 do isym=1,nsym
   call mati3inv(symrel(:,:,isym),symrec(:,:,isym))
 end do

 ! SYMATM generates for all the atoms and all the symmetries, the atom
 ! on which the referenced one is sent and also the translation bringing
 ! back this atom to the referenced unit cell
 ! GA: symatm was already called in crystal_init, no need to do it again.
 call symatm(indsym,natom,nsym,symrec,tnons,tolsym8,typat,xred)

 !write(msg, '(3a,i0,a)' )ch10,ch10,' rdddb9: read ',ddb%nblok,' blocks from the input DDB '
 !call wrtout(std_out,msg)

 ! Read the blocks from the input database, and close it.
 do iblok=1,ddb%nblok

   call ddb%read_block_txt(iblok,mband,mpert,msize,nkpt,unddb)

   if (raw_ == 0) then
     call ddb%symmetrize_and_transform(ddb_hdr%crystal,iblok)
   end if

 end do ! iblok

 DBG_EXIT("COLL")

end subroutine rdddb9
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/chkin9
!! NAME
!! chkin9
!!
!! FUNCTION
!! Check the value of some input parameters.
!! Send error message and stop if needed.
!! Also transform the meaning of atifc
!!
!! INPUTS
!! atifc(natom)=list of the atom ifc to be analysed
!! natifc= number of atom ifc to be analysed
!! natom= number of atoms
!!
!! OUTPUT
!! atifc(natom) =  atifc(ia) equals 1 if the analysis of ifc
!!  has to be done for atom ia; otherwise 0.
!!
!! NOTES
!! Only for one processor (no use of wrtout)
!!
!! SOURCE

subroutine chkin9(atifc,natifc,natom)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: natifc,natom
!arrays
 integer,intent(inout) :: atifc(natom)

!Local variables -------------------------
!scalars
 integer :: iatifc
 character(len=500) :: msg
!arrays
 integer,allocatable :: work(:)

! *********************************************************************

 if(natifc>natom)then
   write(msg, '(a,i0,3a,i0,3a)' )&
    'The number of atom ifc in the input files',natifc,',',ch10,&
    'is larger than the number of atoms',natom,'.',ch10,&
    'Action: change natifc in the input file.'
   ABI_ERROR(msg)
 end if

 if(natifc>=1)then
   ABI_MALLOC(work,(natom))
   work(:)=0

   do iatifc=1,natifc
     if(atifc(iatifc)<=0.or.atifc(iatifc)>natom)then
       write(msg, '(a,i0,5a,i0,3a)' )&
        'For iatifc=',iatifc,', the number of the atom ifc to be ',ch10,&
        'analysed is not valid : either negative, ',ch10,&
        'zero, or larger than natom =',natom,'.',ch10,&
        'Action: change atifc in your input file.'
       ABI_ERROR(msg)
     end if
     work(atifc(iatifc))=1
   end do

   atifc(1:natom)=work(:)
   ABI_FREE(work)
 end if

end subroutine chkin9
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/nlopt
!! NAME
!! nlopt
!!
!! FUNCTION
!! Output of all quantities related to third-order derivatives of the energy.
!! Compute the permutations of the three perturbations, then
!! write out the whole matrix of third order derivatives
!! in reduced coordinates. Finally, compute the non-linear optical
!! susceptibility d and the first-order change in the dielectric
!! susceptibility tensor induced by an atomic displacement.
!!
!! INPUTS
!!  blkflg(3,mpert,3,mpert,3,mpert)= ( 1 if the element of the 3dte
!!   has been calculated ; 0 otherwise )
!!  d3(2,3,mpert,3,mpert,3,mpert)= matrix of the 3DTE
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space(bohr^-1)
!!  mpert =maximum number of ipert
!!  natom= number of atoms
!!  rprimd(3,3)=dimensional primitive translations (bohr)
!!  ucvol=unit cell volume (bohr^3)
!!
!! OUTPUT
!! carflg(3,mpert,3,mpert,3,mpert)=1 if the element of d3cart has been calculated, 0 otherwise
!! d3cart(2,3,mpert,3,mpert,3,mpert)=matrix of third-order energy derivatives in cartesian coordinates
!!
!! SOURCE

subroutine nlopt(blkflg,carflg,d3,d3cart,gprimd,mpert,natom,rprimd,ucvol)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: mpert,natom
 real(dp),intent(in) :: ucvol
!arrays
 integer,intent(in) :: blkflg(3,mpert,3,mpert,3,mpert)
 integer,intent(out) :: carflg(3,mpert,3,mpert,3,mpert)
 real(dp),intent(in) :: d3(2,3,mpert,3,mpert,3,mpert),gprimd(3,3),rprimd(3,3)
 real(dp),intent(out) :: d3cart(2,3,mpert,3,mpert,3,mpert)

!Local variables -------------------------
!scalars
 integer :: i1dir,i1pert,i2dir,i2pert,i3dir,i3pert
!arrays
 integer :: flg1(3),flg2(3)
 real(dp) :: vec1(3),vec2(3)

! *******************************************************************

!Compute the permutations of the perturbations

 d3cart(:,:,:,:,:,:,:) = 0._dp

 do i1pert = 1,mpert
   do i2pert = 1,mpert
     do i3pert = 1,mpert
       do i1dir=1,3
         do i2dir=1,3
           do i3dir=1,3

!            Check if all elements are available

             if ((blkflg(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)/=0).and. &
                 (blkflg(i1dir,i1pert,i3dir,i3pert,i2dir,i2pert)/=0).and. &
                 (blkflg(i2dir,i2pert,i1dir,i1pert,i3dir,i3pert)/=0).and. &
                 (blkflg(i2dir,i2pert,i3dir,i3pert,i1dir,i1pert)/=0).and. &
                 (blkflg(i3dir,i3pert,i1dir,i1pert,i2dir,i2pert)/=0).and. &
                 (blkflg(i3dir,i3pert,i2dir,i2pert,i1dir,i1pert)/=0)) then

               d3cart(:,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = &
               (  d3(:,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) + &
                  d3(:,i1dir,i1pert,i3dir,i3pert,i2dir,i2pert) + &
                  d3(:,i2dir,i2pert,i1dir,i1pert,i3dir,i3pert) + &
                  d3(:,i2dir,i2pert,i3dir,i3pert,i1dir,i1pert) + &
                  d3(:,i3dir,i3pert,i1dir,i1pert,i2dir,i2pert) + &
                  d3(:,i3dir,i3pert,i2dir,i2pert,i1dir,i1pert))*sixth

             end if
           end do
         end do
       end do
     end do
   end do
 end do

!Transform to cartesian coordinates
 carflg(:,:,:,:,:,:) = 0

 do i1pert = 1, mpert
   do i2pert = 1, mpert
     do i3pert = 1, mpert

       do i2dir = 1, 3
         do i3dir = 1, 3

           vec1(:) = d3cart(1,:,i1pert,i2dir,i2pert,i3dir,i3pert)
           flg1(:) = blkflg(:,i1pert,i2dir,i2pert,i3dir,i3pert)
           call cart39(flg1,flg2,gprimd,i1pert,natom,rprimd,vec1,vec2)
           d3cart(1,:,i1pert,i2dir,i2pert,i3dir,i3pert) = vec2(:)
           carflg(:,i1pert,i2dir,i2pert,i3dir,i3pert) = flg2(:)

         end do
       end do

       do i1dir = 1, 3
         do i3dir = 1, 3
           vec1(:) = d3cart(1,i1dir,i1pert,:,i2pert,i3dir,i3pert)
           flg1(:) = blkflg(i1dir,i1pert,:,i2pert,i3dir,i3pert)
           call cart39(flg1,flg2,gprimd,i2pert,natom,rprimd,vec1,vec2)
           d3cart(1,i1dir,i1pert,:,i2pert,i3dir,i3pert) = vec2(:)
           carflg(i1dir,i1pert,:,i2pert,i3dir,i3pert) = flg2(:)
         end do
       end do

       do i1dir = 1, 3
         do i2dir = 1, 3
           vec1(:) = d3cart(1,i1dir,i1pert,i2dir,i2pert,:,i3pert)
           flg1(:) = blkflg(i1dir,i1pert,i2dir,i2pert,:,i3pert)
           call cart39(flg1,flg2,gprimd,i3pert,natom,rprimd,vec1,vec2)
           d3cart(1,i1dir,i1pert,i2dir,i2pert,:,i3pert) = vec2(:)
           carflg(i1dir,i1pert,i2dir,i2pert,:,i3pert) = flg2(:)
         end do
       end do

     end do
   end do
 end do

 ! Compute non linear-optical coefficients d_ijk (atomic units)
 i1pert = natom+2
 d3cart(:,:,i1pert,:,i1pert,:,i1pert) = -3._dp*d3cart(:,:,i1pert,:,i1pert,:,i1pert)/(ucvol*2._dp)

 ! Compute first-order change in the electronic dielectric
 ! susceptibility (Bohr^-1) induced by an atomic displacement
 d3cart(1:2,1:3,1:natom,1:3,natom + 2,1:3,natom + 2) = -6._dp*d3cart(1:2,1:3,1:natom,1:3,natom + 2,1:3,natom + 2)/ucvol

end subroutine nlopt
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_from_file
!! NAME
!!  ddb_from_file
!!
!! FUNCTION
!!  This subroutine reads data from the DDB file and constructs an instance of ddb_type
!!  It also returns an instance of crystal_t with the crystalline structure reported in the DDB file
!!  and the DDB header.
!!
!! INPUTS
!!  filename=DDB filename.
!!  comm=MPI communicator.
!!  [prtvol] = Verbosity level
!!  raw = 1 -> do not perform any symetrization or transformation to cartesian coordinates.
!!        0 (default) -> do perform these transformations.
!!
!! OUTPUT
!!  ddb<type(ddb_type)>=Object storing the DDB results.
!!  crystal<type(crystal_t)>=Crystal structure parameters
!!  ddb_hdr<type(ddb_hdr_type)>= Header of the DDB file.
!!
!! SOURCE

subroutine ddb_from_file(ddb, filename, ddb_hdr, crystal, comm, prtvol, raw)

!Arguments -------------------------------
!scalars
 class(ddb_type),intent(inout) :: ddb
 integer,intent(in) :: comm
 integer,optional,intent(in) :: prtvol, raw
 character(len=*),intent(in) :: filename
 type(crystal_t),intent(out) :: Crystal
 type(ddb_hdr_type),intent(out) :: ddb_hdr
!array

!Local variables-------------------------------
 integer :: iomode
 character(len=fnlen) :: filename_
 integer :: prtvol_
 character(len=500) :: msg

! ************************************************************************

 DBG_ENTER("COLL")

 prtvol_ = 0; if (present(prtvol)) prtvol_ = prtvol

 call ddb_hdr%get_iomode(filename, 1, iomode, filename_)

 if (iomode==IO_MODE_ETSF) then
   call ddb%read_nc(filename_, ddb_hdr, crystal, comm, prtvol, raw)
 else if (iomode==IO_MODE_FORTRAN) then
   call ddb%read_txt(filename_, ddb_hdr, crystal, comm, prtvol, raw)
 end if

 ! Print out info on the crystal
 if (prtvol_ >= 0) then

   call ddb_hdr%crystal%print(unit=ab_out)
   call ddb_hdr%crystal%print(unit=std_out)

   write(msg, '(2a,i0,a)' )ch10,' DDB file with ',ddb%nblok,' blocks has been read.'
   call wrtout(std_out,msg)
   call wrtout(ab_out,msg)

 end if

 DBG_EXIT("COLL")

end subroutine ddb_from_file
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_read_txt
!! NAME
!!  ddb_read_txt
!!
!! FUNCTION
!!  This subroutine reads data from the DDB file and constructs an instance of ddb_type
!!  It also returns an instance of crystal_t with the crystalline structure reported in the DDB file
!!  and the DDB header.
!!
!! INPUTS
!!  filename=DDB filename.
!!  comm=MPI communicator.
!!  prtvol=Verbosity level
!!  raw = 1 -> do not perform any symetrization or transformation to cartesian coordinates.
!!        0 (default) -> do perform these transformations.
!!
!! OUTPUT
!!  ddb<type(ddb_type)>=Object storing the DDB results.
!!  crystal<type(crystal_t)>=Crystal structure parameters
!!  ddb_hdr= Header of the DDB file.
!!
!! SOURCE

subroutine ddb_read_txt(ddb, filename, ddb_hdr, crystal, comm, prtvol, raw)

!Arguments -------------------------------
!scalars
 class(ddb_type),intent(inout) :: ddb
 integer,intent(in) :: comm
 integer,optional,intent(in) :: prtvol, raw
 character(len=*),intent(in) :: filename
 type(crystal_t),intent(out) :: Crystal
 type(ddb_hdr_type),intent(out) :: ddb_hdr

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: msym,dimekb,lmnmax,mband,nkpt,ntypat,nsym,usepaw
 integer :: mpert,msize,natom,nblok,occopt,nsppol
 real(dp) :: ucvol

!arrays
 integer,allocatable :: symrec(:,:,:),symrel(:,:,:),symafm(:),indsym(:,:,:),typat(:)
 real(dp) :: acell(3),gmet(3,3),gprim(3,3),rmet(3,3),rprim(3,3)
 real(dp),allocatable :: amu(:),xcart(:),xred(:,:),zion(:),znucl(:),tnons(:,:)

! ************************************************************************

 DBG_ENTER("COLL")

! Must read natom from the DDB before being able to allocate some arrays needed for invars9
 call ddb_hdr%open_read_txt(filename, comm)

! GA: clean this up. Not all of it is useful
 nblok = ddb_hdr%nblok
 msym = ddb_hdr%msym
 natom = ddb_hdr%natom
 ntypat = ddb_hdr%ntypat
 mband = ddb_hdr%mband
 nkpt = ddb_hdr%nkpt
 nsppol = ddb_hdr%nsppol
 usepaw = ddb_hdr%usepaw
 dimekb = ddb_hdr%psps%dimekb
 lmnmax = ddb_hdr%psps%lmnmax
 occopt = ddb_hdr%occopt
 mpert = ddb_hdr%mpert
 msize = ddb_hdr%msize

 ! Master reads and then broadcasts data.
 if (xmpi_comm_rank(comm) == master) then

   ! Allocate arrays depending on msym (which is actually fixed to nsym inside inprep8)
   ABI_MALLOC(symrel,(3,3,msym))
   ABI_MALLOC(symafm,(msym))
   ABI_MALLOC(tnons,(3,msym))
   ABI_MALLOC(typat,(natom))
   ABI_MALLOC(xred,(3,natom))
   ABI_MALLOC(zion,(ntypat))
   ABI_MALLOC(znucl,(ntypat))

   ABI_MALLOC(symrec,(3,3,msym))
   ABI_MALLOC(indsym,(4,msym,natom))
   ABI_MALLOC(xcart,(3*natom))
   ABI_MALLOC(amu,(ntypat))

   ddb%nsppol = nsppol
   call ddb%malloc(msize, nblok, natom, ntypat, mpert)

   ! GA: FIXME
   ! Should clean this up. Lots of the arguments are not needed.
   ! In particular, rprim and acell could be taken from ddb_hdr%crystal
   ! which is already initialized at this point
   call rdddb9(ddb, ddb_hdr, ddb_hdr%unddb,&
    acell,amu,gmet,gprim,indsym,&
    mband,mpert,msize,msym,&
    natom,nkpt,nsym,ntypat,&
    rmet,rprim,symrec,symrel,symafm,&
    tnons,typat,ucvol,xcart,xred,zion,znucl,raw)

   ABI_FREE(symrec)
   ABI_FREE(indsym)
   ABI_FREE(xcart)

   ! Save variables needed to call legacy code.
   ddb%acell = acell
   ddb%rprim = rprim
   ddb%gprim = gprim

   !call ddb%set_brav(brav)

   ! Other useful quantities.
   ! 2 is to preserve the old behaviour
   ddb%prtvol = 2; if (present(prtvol)) ddb%prtvol = prtvol
   ddb%occopt = occopt
   ddb%amu = amu
   ABI_FREE(amu)

   ! These were not needed because crystal is already initialized in the header.
   ABI_FREE(symrel)
   ABI_FREE(symafm)
   ABI_FREE(tnons)
   ABI_FREE(typat)
   ABI_FREE(xred)
   ABI_FREE(zion)
   ABI_FREE(znucl)

 end if

 call ddb_hdr%close()

 if (xmpi_comm_size(comm) > 1) then
   call ddb%bcast(comm)
   call ddb_hdr%bcast(comm)

   !! GA: This seems superfluous now...
   !call xmpi_bcast(nsym, master, comm, ierr)
   !call xmpi_bcast(symrel, master, comm, ierr)
   !call xmpi_bcast(symafm, master, comm, ierr)
   !call xmpi_bcast(typat, master, comm, ierr)
   !call xmpi_bcast(acell, master, comm, ierr)
   !call xmpi_bcast(occopt, master, comm, ierr)
   !call xmpi_bcast(gprim, master, comm, ierr)
   !call xmpi_bcast(rprim, master, comm, ierr)
   !call xmpi_bcast(tnons, master, comm, ierr)
   !call xmpi_bcast(xred, master, comm, ierr)
   !call xmpi_bcast(zion, master, comm, ierr)
   !call xmpi_bcast(znucl, master, comm, ierr)
 end if

 call ddb_hdr%crystal%copy(Crystal)

 !! Initialize crystal_t object.
 !call mkrdim(acell,rprim,rprimd)

 !! GA: These variables are hardcoded which means the crystal object
 !!     is not reliable for antiferro systems or alchemical potentials
 !!     when it is read from a text DDB file.
 !npsp = ntypat; space_group = 0; timrev = 2
 !use_antiferro=.FALSE. !;  use_antiferro=(nspden==2.and.nsppol==1)
 !ABI_MALLOC(title, (ntypat))

 !do ii=1,ntypat
 !  write(title(ii),'(a,i0)')"No title for typat ",ii
 !end do

 !! Warning znucl is dimensioned with ntypat = nspsp hence alchemy is not supported here
 !call crystal_init(ddb%amu,Crystal,space_group,natom,npsp,ntypat,nsym,rprimd,typat,xred,&
 !  zion,znucl,timrev,use_antiferro,.FALSE.,title,&
 !  symrel=symrel(:,:,1:nsym),tnons=tnons(:,1:nsym),symafm=symafm(1:nsym))

 !ABI_FREE(title)
 !ABI_FREE(symrel)
 !ABI_FREE(symafm)
 !ABI_FREE(tnons)
 !ABI_FREE(typat)
 !ABI_FREE(xred)
 !ABI_FREE(zion)
 !ABI_FREE(znucl)

 DBG_EXIT("COLL")

end subroutine ddb_read_txt
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_read_nc
!! NAME
!!  ddb_read_nc
!!
!! FUNCTION
!!  This subroutine reads data from the DDB.nc file and constructs an instance of ddb_type
!!  It also returns an instance of crystal_t with the crystalline structure reported in the DDB file
!!  and the DDB header.
!!
!! INPUTS
!!  filename=DDB filename.
!!  comm=MPI communicator.
!!  prtvol=Verbosity level
!!  raw = 1 -> do not perform any symetrization or transformation to cartesian coordinates.
!!        0 (default) -> do perform these transformations.
!!
!! OUTPUT
!!  ddb<type(ddb_type)>=Object storing the DDB results.
!!  crystal<type(crystal_t)>=Crystal structure parameters
!!  ddb_hdr= Header of the DDB file.
!!
!! SOURCE

subroutine ddb_read_nc(ddb, filename, ddb_hdr, crystal, comm, prtvol, raw)

!Arguments -------------------------------
!scalars
 class(ddb_type),intent(inout) :: ddb
 integer,intent(in) :: comm
 integer,optional,intent(in) :: prtvol, raw
 character(len=*),intent(in) :: filename
 type(crystal_t),intent(out) :: crystal
 type(ddb_hdr_type),intent(out) :: ddb_hdr
!array

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: prtvol_, raw_
 integer :: ncid
 integer :: iblok,iblok_d0E,iblok_d1E,iblok_d2E,iblok_d3E,iblok_d2eig

!arrays
 !character(len=132),allocatable :: title(:)

! ************************************************************************

 DBG_ENTER("COLL")

 if (present(raw)) then
   raw_ = raw
 else
   raw_ = 0
 end if

 ! GA: Not really used so far
 if (present(prtvol)) then
   prtvol_ = prtvol
 else
   prtvol_ = 0
 end if

 ! Read header
 call ddb_hdr%open_read_nc(filename, comm)
 ncid = ddb_hdr%ncid

 if (xmpi_comm_rank(comm) == master) then

   ddb%nsppol = ddb_hdr%nsppol

   ! Copy dimensions from header and allocate arrays
   call ddb%malloc(ddb_hdr%msize, ddb_hdr%nblok, ddb_hdr%natom, &
                   ddb_hdr%ntypat, ddb_hdr%mpert,&
                   ddb_hdr%nkpt, ddb_hdr%mband)

   ! Copy arrays from header
   ddb%typ(:) = ddb_hdr%typ(:)
   ddb%amu(:) = ddb_hdr%crystal%amu(:)
   ddb%acell(:) = one
   ddb%rprim(:,:) = ddb_hdr%crystal%rprimd(:,:)
   ddb%gprim(:,:) = ddb_hdr%crystal%gprimd(:,:)

   ! ---------------
   ! Read all blocks
   ! ---------------
   iblok_d0E = 0
   iblok_d1E = 0
   iblok_d2E = 0
   iblok_d3E = 0
   iblok_d2eig = 0

   do iblok=1,ddb%nblok

     if (is_type_d0E(ddb%typ(iblok))) then
       iblok_d0E = iblok_d0E + 1
       call ddb%read_d0E_nc(ncid, iblok, iblok_d0E)

     else if (is_type_d1E(ddb%typ(iblok))) then
       iblok_d1E = iblok_d1E + 1
       call ddb%read_d1E_nc(ncid, iblok, iblok_d1E)

     else if (is_type_d2E(ddb%typ(iblok))) then
       iblok_d2E = iblok_d2E + 1
       call ddb%read_d2E_nc(ncid, iblok, iblok_d2E)

     else if (is_type_d3E(ddb%typ(iblok))) then
       iblok_d3E = iblok_d3E + 1
       call ddb%read_d3E_nc(ncid, iblok, iblok_d3E)

     else if (is_type_d2eig(ddb%typ(iblok))) then
       iblok_d2eig = iblok_d2eig + 1
       ! GA: It is kind of weird to call this function inside a loop,
       !     because the ddb can only hold a single block of d2eig data.
       call ddb%read_d2eig_nc(ncid, iblok, iblok_d2eig)

     end if

     ! Symmetrize and transform if raw==0
     if (raw_==0) then
       call ddb%symmetrize_and_transform(ddb_hdr%crystal,iblok)
     end if

   end do

 end if

 ! Close the file
 call ddb_hdr%close()

 ! --------------
 ! Broadcast data
 ! --------------
 if (xmpi_comm_size(comm) > 1) then
   call ddb%bcast(comm)
   call ddb_hdr%bcast(comm)
 end if

 ! Copy crystal
 call ddb_hdr%crystal%copy(crystal)

 DBG_EXIT("COLL")

end subroutine ddb_read_nc
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_can_merge_blocks
!! NAME
!! ddb_can_merge_blocks
!!
!! FUNCTION
!!  Return true if iblok1 of ddb1 can be merged to iblok2 of ddb2
!!
!! INPUTS
!!  ddb1=ddb object 1
!!  ddb2=ddb object 2
!!  iblok1=block index from ddb1
!!  iblok2=block index from ddb2
!!
!! OUTPUT
!!  can_merge=.true. if the blocks are compatible for merging.
!!
!! SOURCE

logical function ddb_can_merge_blocks(ddb1, ddb2, iblok1, iblok2) result(can_merge)

!Arguments -------------------------------
!array
 class(ddb_type),intent(inout) :: ddb1
 type(ddb_type),intent(inout) :: ddb2
 integer,intent(in) :: iblok1
 integer,intent(in) :: iblok2

!local variables
!scalars
 integer :: nq, ii, blktyp
 real(dp),parameter :: qtol=2.0d-8
 real(dp) :: diff

! ************************************************************************

  can_merge = .false.
  if(ddb1%typ(iblok1)/=ddb2%typ(iblok2)) return

  blktyp = ddb1%typ(iblok1)

  can_merge = .true.

  if (is_type_d0E(blktyp) .or. is_type_d1E(blktyp)) return

  ! Compare wavevectors
  if (is_type_d2E(blktyp).or.is_type_d2eig(blktyp))then
    nq=1
  else if (is_type_d3E(blktyp))then
    nq=3
  end if

  do ii=1,nq
    diff = (ddb1%qpt(1+3*(ii-1),iblok1)/ddb1%nrm(ii,iblok1) &
          - ddb2%qpt(1+3*(ii-1),iblok2)/ddb2%nrm(ii,iblok2))
    if (abs(diff) > qtol) can_merge = .false.
    diff = (ddb1%qpt(2+3*(ii-1),iblok1)/ddb1%nrm(ii,iblok1) &
          - ddb2%qpt(2+3*(ii-1),iblok2)/ddb2%nrm(ii,iblok2))
    if (abs(diff) > qtol) can_merge = .false.
    diff = (ddb1%qpt(3+3*(ii-1),iblok1)/ddb1%nrm(ii,iblok1) &
          - ddb2%qpt(3+3*(ii-1),iblok2)/ddb2%nrm(ii,iblok2))
    if (abs(diff) > qtol) can_merge = .false.
  end do

end function ddb_can_merge_blocks
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_merge_blocks
!! NAME
!! ddb_merge_blocks
!!
!! FUNCTION
!!  Merge block number iblok2 from ddb2 into block number iblok1 in ddb1.
!!
!! INPUTS
!!  ddb1=ddb object 1
!!  ddb2=ddb object 2
!!  iblok1=block index from ddb1
!!  iblok2=block index from ddb2
!!
!! OUTPUT
!!
!! SOURCE

subroutine ddb_merge_blocks(ddb1, ddb2, iblok1, iblok2)

!Arguments -------------------------------
!array
 class(ddb_type),intent(inout) :: ddb1
 type(ddb_type),intent(inout) :: ddb2
 integer,intent(in) :: iblok1
 integer,intent(in) :: iblok2

!local variables
!scalars
 integer :: ii, blktyp, mpert1, mpert2
 integer :: idir1, idir2, idir3, ipert1, ipert2, ipert3
 real(dp),parameter :: qtol=2.0d-8
!arrays
 real(dp), allocatable :: d1matr(:,:,:)
 real(dp), allocatable :: d2matr(:,:,:,:,:)
 real(dp), allocatable :: d3matr(:,:,:,:,:,:,:)
 integer, allocatable  :: d1flg(:,:)
 integer, allocatable  :: d2flg(:,:,:,:)
 integer, allocatable  :: d3flg(:,:,:,:,:,:)

! ************************************************************************

  ! Note that ddb and ddb2 may have a different values for mpert
  !mpert = min(ddb1%mpert, ddb2%mpert)
  mpert1 = ddb1%mpert
  mpert2 = ddb2%mpert

  ! Add the blok to the output ddb
  blktyp = ddb2%typ(iblok2)
  ddb1%typ(iblok1) = blktyp

  ! Copy q-point
  do ii=1,9
    ddb1%qpt(ii,iblok1) = ddb2%qpt(ii,iblok2)
  end do
  do ii=1,3
    ddb1%nrm(ii,iblok1) = ddb2%nrm(ii,iblok2)
  end do
  
  if (is_type_d0E(blktyp)) then
     ! --------------
     ! Copy d0E block
     ! --------------
     if (ddb2%flg(1,iblok2) > 0) then
       ddb1%val(1,1,iblok1) = ddb2%val(1,1,iblok2)
       ddb1%val(2,1,iblok1) = ddb2%val(2,1,iblok2)
       ddb1%flg(1,iblok1) = ddb2%flg(1,iblok2)
     end if

  else if (is_type_d1E(blktyp)) then
     ! --------------
     ! Copy d1E block
     ! --------------
     call ddb2%get_d1matr(iblok2, d1matr, d1flg)
     !call ddb%set_d1matr(iblok, d1matr, d1flg)
     ii=0
     do ipert1=1,mpert1
       do idir1=1,3
         ii=ii+1
         if (ipert1 <= mpert2) then
           if (d1flg(idir1,ipert1)>0) then
             ddb1%val(1,ii,iblok1) = d1matr(1,idir1,ipert1)
             ddb1%val(2,ii,iblok1) = d1matr(2,idir1,ipert1)
             ddb1%flg(ii,iblok1) = d1flg(idir1,ipert1)
           end if
         end if
       end do
     end do
     ABI_SFREE(d1matr)
     ABI_SFREE(d1flg)

  else if (is_type_d2E(blktyp)) then
     ! --------------
     ! Copy d2E block
     ! --------------
     call ddb2%get_d2matr(iblok2, d2matr, d2flg)
     !call ddb%set_d2matr(iblok, d2matr, d2flg)
     ii=0
     do ipert2=1,mpert1
       do idir2=1,3
         do ipert1=1,mpert1
           do idir1=1,3
             ii=ii+1
             if ((ipert1 <= mpert2).and.(ipert2<=mpert2)) then
               if (d2flg(idir1,ipert1,idir2,ipert2)>0) then
                 ddb1%val(1,ii,iblok1) = d2matr(1,idir1,ipert1,idir2,ipert2)
                 ddb1%val(2,ii,iblok1) = d2matr(2,idir1,ipert1,idir2,ipert2)
                 ddb1%flg(ii,iblok1) = d2flg(idir1,ipert1,idir2,ipert2)
               end if
             end if
           end do
         end do
       end do
     end do
     ABI_SFREE(d2matr)
     ABI_SFREE(d2flg)

  else if (is_type_d3E(blktyp)) then
     ! --------------
     ! Copy d3E block
     ! --------------
     call ddb2%get_d3matr(iblok2, d3matr, d3flg)
     !call ddb%set_d3matr(iblok, d3matr, d3flg)
     ii=0
     do ipert1=1,mpert1
       do idir1=1,3
         do ipert2=1,mpert1
           do idir2=1,3
             do ipert3=1,mpert1
               do idir3=1,3


                 !ii=ii+1  ! GA: This is not equivalent
                 ii = idir1 + 3*((ipert1-1)+mpert1*((idir2-1) &
                            + 3*((ipert2-1)+mpert1*((idir3-1) &
                            + 3*(ipert3-1)))))

                 ! Note that the loop order (1,2,3) is not really consistent
                 ! with the d2E case (2, 1)
                 ! TODO Clean this up

                 if ((ipert1 <= mpert2).and.(ipert2<=mpert2).and.(ipert3<=mpert2)) then
                   if (d3flg(idir1,ipert1,idir2,ipert2,idir3,ipert3)>0) then
                     ddb1%flg(ii,iblok1) = d3flg(idir1,ipert1,idir2,ipert2,idir3,ipert3)
                     ddb1%val(:,ii,iblok1) = d3matr(:,idir1,ipert1,idir2,ipert2,idir3,ipert3)
                   end if
                 end if

               end do
             end do
           end do
         end do
       end do
     end do
     ABI_SFREE(d3matr)
     ABI_SFREE(d3flg)

  else if (is_type_d2eig(blktyp)) then
     ! ----------------
     ! Skip d2eig block
     ! ----------------

     ! TODO need a function ddb_merge_d2eig(filename1, filename2, )

  end if  ! blktyp

end subroutine ddb_merge_blocks
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/carttransf
!! NAME
!! carttransf
!!
!! FUNCTION
!! Transform a second-derivative matrix (EIG2D) from reduced
!! coordinates to cartesian coordinates.
!!
!! INPUTS
!!  blkflg(msize,nblok)=
!!   ( 1 if the element of the dynamical matrix has been calculated ;
!!     0 otherwise )
!!  gprimd(3,3)=basis vector in the reciprocal space
!!  iqpt  = number of the Q-point currently used
!!  mband = maximal number of bands
!!  mpert = maximum number of ipert
!!  msize = size of the EIG2D arrays (3*mpert*3*mpert)
!!  natom = number of atom
!!  nblok = number of bloks in blkflg
!!  nkpt  = number of K-points
!!  rprimd(3,3) = basis vector in the real space
!!
!! OUTPUT
!!  carflg(3,mpert,3,mpert)= ( 1 if the element of the cartesian
!!  EIG2D matrix has been calculated correctly ; 0 otherwise )
!!
!! SIDE EFFECT
!! blkval2(2,msize,mband,nkpt)=Second order eigenvalues (EIG2D)
!! is transformed from reduced coordinates to cartesian coordinates
!!
!! SOURCE

subroutine carttransf(blkflg,blkval2,carflg,gprimd,iqpt,mband, mpert,msize,natom,nblok,nkpt,rprimd)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: mband,msize
 integer,intent(in) :: iqpt
 integer,intent(in) :: mpert,nblok
 integer,intent(inout) :: natom,nkpt
!arrays
 integer,intent(in) :: blkflg(msize,nblok)
 integer,intent(out) :: carflg(3,mpert,3,mpert)
 real(dp),intent(in) :: gprimd(3,3),rprimd(3,3)
 real(dp),intent(inout) :: blkval2(2,msize,mband,nkpt)

!Local variables-------------------------------
!scalars
integer :: iatom1,iatom2,iband,idir1,idir2,ikpt, index
!arrays
real(dp),allocatable :: blkflgtmp(:,:,:,:,:), blkval2tmp(:,:,:,:,:,:), d2cart(:,:,:,:,:)

! *********************************************************************

 ! Start by allocating local arrays
 ABI_MALLOC(blkflgtmp,(3,mpert,3,mpert,1))
 ABI_MALLOC(blkval2tmp,(2,3,mpert,3,mpert,1))
 ABI_MALLOC(d2cart,(2,3,mpert,3,mpert))

 ! Begin by formating the arrays to be compatible with cart29
 ! Then call cart29 to transform the arrays in cartesian coordinates
 ! Finally reformat the cartesian arrays in old format
 do ikpt=1,nkpt
   do iband=1,mband

     do idir1=1,3
       do iatom1=1,mpert
         do idir2=1,3
           do iatom2=1,mpert
             index = idir1 + 3*((iatom1 - 1) + natom * ((idir2-1)+3*(iatom2-1)))
             blkflgtmp(idir1,iatom1,idir2,iatom2,1) = blkflg(index,iqpt)
             blkval2tmp(:,idir1,iatom1,idir2,iatom2,1) = blkval2(:,index,iband,ikpt)
           end do
         end do
       end do
     end do

     ! The 1sin the argument of cart29 are respectively iblok and nblok. We are doing only one blok.
     call carteig2d(blkflg(:,iqpt),blkval2tmp,carflg,d2cart,gprimd,1,mpert,natom,1,rprimd)

     do idir1=1,3
       do iatom1=1,mpert
         do idir2=1,3
           do iatom2=1,mpert
             index = idir1 + 3*((iatom1 - 1) + natom * ((idir2-1)+3*(iatom2-1)))
             blkval2(:,index,iband,ikpt) = d2cart(:,idir1,iatom1,idir2,iatom2)
           end do
         end do
       end do
     end do

   end do
 end do

 ABI_FREE(blkflgtmp)
 ABI_FREE(blkval2tmp)
 ABI_FREE(d2cart)

end subroutine carttransf
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/carteig2d
!! NAME
!! carteig2d
!!
!! FUNCTION
!! Transform a second-derivative matrix (EIG2D) from reduced
!! coordinates to cartesian coordinates
!!
!! INPUTS
!!  blkflg(3,mpert,3,mpert,nblok)=
!!   ( 1 if the element of the dynamical matrix has been calculated ;
!!     0 otherwise )
!!  blkval(2,3,mpert,3,mpert,nblok)=DDB values
!!  gprimd(3,3)=basis vector in the reciprocal space
!!  iblok=number of the blok that will be transformed
!!  mpert =maximum number of ipert
!!  natom=number of atom
!!  nblok=number of blocks (dimension of blkflg and blkval)
!!  rprimd(3,3)=basis vector in the real space
!!
!! OUTPUT
!!  carflg(3,mpert,3,mpert)= ( 1 if the element of the cartesian
!!  2DTE matrix has been calculated correctly ; 0 otherwise )
!!  d2cart(2,3,mpert,3,mpert)=
!!    dynamical matrix, effective charges, dielectric tensor,....
!!    all in cartesian coordinates
!!
!! SOURCE

subroutine carteig2d(blkflg,blkval,carflg,d2cart,gprimd,iblok,mpert,natom,nblok,rprimd)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iblok,mpert,natom,nblok
!arrays
 integer,intent(in) :: blkflg(3,mpert,3,mpert,nblok)
 integer,intent(out) :: carflg(3,mpert,3,mpert)
 real(dp),intent(in) :: blkval(2,3,mpert,3,mpert,nblok),gprimd(3,3),rprimd(3,3)
 real(dp),intent(out) :: d2cart(2,3,mpert,3,mpert)

!Local variables -------------------------
!scalars
 integer :: idir1,idir2,ii,ipert1,ipert2
!arrays
 integer :: flg1(3),flg2(3)
 real(dp) :: vec1(3),vec2(3)

! *********************************************************************

 ! First, copy the data blok in place.
 d2cart(:,:,:,:,:)=blkval(:,:,:,:,:,iblok)

 ! Cartesian coordinates transformation (in two steps)
 ! First step
 do ipert1=1,mpert
   do ipert2=1,mpert
     do ii=1,2
       do idir1=1,3
         do idir2=1,3
           vec1(idir2)=d2cart(ii,idir1,ipert1,idir2,ipert2)
           ! Note here blkflg
           flg1(idir2)=blkflg(idir1,ipert1,idir2,ipert2,iblok)
         end do
         call cart39(flg1,flg2,gprimd,ipert2,natom,rprimd,vec1,vec2)
         do idir2=1,3
           d2cart(ii,idir1,ipert1,idir2,ipert2)=vec2(idir2)
           ! And here carflg
           carflg(idir1,ipert1,idir2,ipert2)=flg2(idir2)
         end do
       end do
     end do
   end do
 end do

 ! Second step
 do ipert1=1,mpert
   do ipert2=1,mpert
     do ii=1,2
       do idir2=1,3
         do idir1=1,3
           vec1(idir1)=d2cart(ii,idir1,ipert1,idir2,ipert2)
           ! Note here carflg
           flg1(idir1)=carflg(idir1,ipert1,idir2,ipert2)
         end do
         call cart39(flg1,flg2,gprimd,ipert1,natom,rprimd,vec1,vec2)
         do idir1=1,3
           d2cart(ii,idir1,ipert1,idir2,ipert2)=vec2(idir1)
           ! And here carflg again
           carflg(idir1,ipert1,idir2,ipert2)=flg2(idir1)
         end do
       end do
     end do
   end do
 end do

end subroutine carteig2d
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/dtech9
!!
!! NAME
!! dtech9
!!
!! FUNCTION
!! Reads the Dielectric Tensor and the Effective Charges in the
!! Gamma Block coming from the Derivative Data Base.
!!
!! INPUTS
!! natom= number of atoms in unit cell
!! iblok= index of the Gamma block
!! mpert =maximum number of ipert
!! nblok= number of blocks in the DDB
!! blkval(2,3*mpert*3*mpert,nblok)=  dynamical matrices
!!  In our case, the nblok is restricted to iblok
!! [unit]=Output unit number
!!
!! OUTPUT
!! zeff(3,3,natom)=effective charge on each atom, versus electric
!!  field and atomic displacement. Note the following convention:
!!  zeff(electric field direction, atomic direction, atom index)
!! dielt(3,3)=dielectric tensor
!!
!! SOURCE

subroutine dtech9(blkval,dielt,iblok,mpert,natom,nblok,zeff,unit)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iblok,mpert,natom,nblok
 integer,intent(in),optional :: unit
!arrays
 real(dp),intent(in) :: blkval(2,3,mpert,3,mpert,nblok)
 real(dp),intent(out) :: dielt(3,3),zeff(3,3,natom)

!Local variables -------------------------
!scalars
 integer :: depl,elec,elec1,elec2,iatom, unt
 character(len=1000) :: msg

! *********************************************************************

 unt = std_out; if (present(unit)) unt = unit

 ! Extraction of effectives charges
 do iatom=1,natom
   do elec=1,3
     do depl=1,3
       zeff(elec,depl,iatom)=0.5*&
        (blkval(1,depl,iatom,elec,natom+2,iblok)+&
         blkval(1,elec,natom+2,depl,iatom,iblok))
     end do
   end do
 end do

 ! Extraction of dielectric tensor
 do elec1=1,3
   do elec2=1,3
     dielt(elec1,elec2)=blkval(1,elec1,natom+2,elec2,natom+2,iblok)
   end do
 end do

 write(msg,'(a,3es16.6,3es16.6,3es16.6)' )' Dielectric Tensor ',&
   dielt(1,1),dielt(1,2),dielt(1,3),&
   dielt(2,1),dielt(2,2),dielt(2,3),&
   dielt(3,1),dielt(3,2),dielt(3,3)
 call wrtout(unt, msg)

 call wrtout(unt, ' Effectives Charges ')
 do iatom=1,natom
   write(msg,'(a,i4,3es16.6,3es16.6,3es16.6)' )' atom ',iatom,&
     zeff(1,1,iatom),zeff(1,2,iatom),zeff(1,3,iatom),&
     zeff(2,1,iatom),zeff(2,2,iatom),zeff(2,3,iatom),&
     zeff(3,1,iatom),zeff(3,2,iatom),zeff(3,3,iatom)
    call wrtout(unt, msg)
 end do

end subroutine dtech9
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/dtchi
!! NAME
!! dtchi
!!
!! FUNCTION
!! Reads the non-linear optical susceptibility tensor and the
!! first-order change in the linear dielectric susceptibility
!! induced by an atomic displacement in the Gamma Block coming from the Derivative Data Base
!! (third-order derivatives).
!!
!! INPUTS
!! blkval(2,3*mpert*3*mpert*3*mpert)= matrix of third-order energies
!! natom= number of atoms in unit cell
!! mpert =maximum number of ipert
!! ramansr= if /= 0, impose sum rule on first-order derivatives
!!                   of the electronic susceptibility with respect
!!                   to atomic displacements
!! nlflag= if =3, only the non-linear optical susceptibilities is computed
!!
!! OUTPUT
!! dchide(3,3,3) = non-linear optical coefficients
!! dchidt(natom,3,3,3) = first-order change of the electronic dielectric
!!   tensor induced by an individual atomic displacement
!!
!! SOURCE

subroutine dtchi(blkval,dchide,dchidt,mpert,natom,ramansr,nlflag)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: mpert,natom,ramansr,nlflag
!arrays
 real(dp),intent(in) :: blkval(2,3*mpert*3*mpert*3*mpert)
 real(dp),intent(out) :: dchide(3,3,3),dchidt(natom,3,3,3)

!Local variables -------------------------
!scalars
 integer :: depl,elfd1,elfd2,elfd3,iatom,ivoigt
 logical :: iwrite
 real(dp) :: wttot
!arrays
 integer :: voigtindex(6,2)
 real(dp) :: d3cart(2,3,mpert,3,mpert,3,mpert),dvoigt(3,6),sumrule(3,3,3)
 real(dp) :: wghtat(natom)

! *********************************************************************

 d3cart(1,:,:,:,:,:,:) = reshape(blkval(1,:),shape = (/3,mpert,3,mpert,3,mpert/))
 d3cart(2,:,:,:,:,:,:) = reshape(blkval(2,:),shape = (/3,mpert,3,mpert,3,mpert/))

 ! Extraction of non-linear optical coefficients
 do elfd1 = 1,3
   do elfd2 = 1,3
     do elfd3 = 1,3
       dchide(elfd1,elfd2,elfd3) = d3cart(1,elfd1,natom+2,elfd2,natom+2,elfd3,natom+2)
     end do
   end do
 end do

 ! Transform to Voigt notations
 voigtindex(:,1) = (/1,2,3,2,1,1/)
 voigtindex(:,2) = (/1,2,3,3,3,2/)
 do ivoigt = 1, 6
   elfd2 = voigtindex(ivoigt,1)
   elfd3 = voigtindex(ivoigt,2)
   do elfd1 = 1, 3
     dvoigt(elfd1,ivoigt) = 0.5_dp*(dchide(elfd1,elfd2,elfd3) + dchide(elfd1,elfd3,elfd2))
   end do
 end do

 ! Transform to pm/V
 dvoigt(:,:) = dvoigt(:,:)*16*(pi**2)*(Bohr_Ang**2)*1.0d-8*eps0/e_Cb

 ! Extraction of $\frac{d \chi}{d \tau}$
 if (nlflag < 3) then
   do iatom = 1, natom
     do depl = 1,3
       do elfd1 = 1,3
         do elfd2 = 1,3
           dchidt(iatom,depl,elfd1,elfd2) = d3cart(1,depl,iatom,elfd1,natom+2,elfd2,natom+2)
         end do
       end do
     end do
   end do
 end if

 wghtat(:) = zero
 if (ramansr == 1) then
   wghtat(:) = one/dble(natom)

 else if (ramansr == 2) then

   wttot = zero
   do iatom = 1, natom
     do depl = 1,3
       do elfd1 = 1,3
         do elfd2 = 1,3
           wghtat(iatom) = wghtat(iatom) + abs(dchidt(iatom,depl,elfd1,elfd2))
         end do
       end do
     end do
     wttot = wttot + wghtat(iatom)
   end do

   wghtat(:) = wghtat(:)/wttot
 end if

 iwrite = ab_out > 0

 if (iwrite) then
   write(ab_out,*)ch10
   write(ab_out,*)'Non-linear optical coefficients d (pm/V)'
   write(ab_out,'(6f12.6)')dvoigt(1,:)
   write(ab_out,'(6f12.6)')dvoigt(2,:)
   write(ab_out,'(6f12.6)')dvoigt(3,:)
 end if

 if (ramansr /= 0) then
   if (iwrite) then
     write(ab_out,*)ch10
     write(ab_out,*)'The violation of the Raman sum rule'
     write(ab_out,*)'by the first-order electronic dielectric tensors ','is as follows'
     write(ab_out,*)'    atom'
     write(ab_out,*)' displacement'
   end if

   sumrule(:,:,:) = zero
   do elfd2 = 1,3
     do elfd1 = 1,3
       do depl = 1,3
         do iatom = 1, natom
           sumrule(depl,elfd1,elfd2) = sumrule(depl,elfd1,elfd2) + dchidt(iatom,depl,elfd1,elfd2)
         end do
         do iatom = 1, natom
           dchidt(iatom,depl,elfd1,elfd2) = dchidt(iatom,depl,elfd1,elfd2) - wghtat(iatom)*sumrule(depl,elfd1,elfd2)
         end do
       end do
     end do
   end do

   if (iwrite) then
     do depl = 1,3
       write(ab_out,'(6x,i2,3(3x,f16.9))') depl,sumrule(depl,1,1:3)
       write(ab_out,'(8x,3(3x,f16.9))') sumrule(depl,2,1:3)
       write(ab_out,'(8x,3(3x,f16.9))') sumrule(depl,3,1:3)
       write(ab_out,*)
     end do
    end if
 end if    ! ramansr

 if (nlflag < 3) then
   if (iwrite) then
     write(ab_out,*)ch10
     write(ab_out,*)' First-order change in the electronic dielectric '
     write(ab_out,*)' susceptibility tensor (Bohr^-1)'
     write(ab_out,*)' induced by an atomic displacement'
     if (ramansr /= 0) then
       write(ab_out,*)' (after imposing the sum over all atoms to vanish)'
     end if
     write(ab_out,*)'  atom  displacement'

     do iatom = 1,natom
       do depl = 1,3
         write(ab_out,'(1x,i4,9x,i2,3(3x,f16.9))')iatom,depl,dchidt(iatom,depl,1,:)
         write(ab_out,'(16x,3(3x,f16.9))')dchidt(iatom,depl,2,:)
         write(ab_out,'(16x,3(3x,f16.9))')dchidt(iatom,depl,3,:)
       end do

       write(ab_out,*)
     end do
   end if
 end if

end subroutine dtchi
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_get_etotal
!!
!! NAME
!!  ddb_get_etotal
!!
!! FUNCTION
!!  Read the GS total energy from the DDB file.
!!
!! INPUTS
!!  ddb<type(ddb_type)>=Derivative database.
!!
!! OUTPUT
!!  etotal=GS Total energy in Hartree
!!  iblok=Index of the block in the DDB file. 0 if not found.
!!
!! SOURCE

integer function ddb_get_etotal(ddb, etotal) result(iblok)

!Arguments -------------------------------
!scalars
 real(dp),intent(out) :: etotal
 class(ddb_type),intent(in) :: ddb

!Local variables -------------------------
!scalars
 integer :: rftyp
!arrays
 integer :: rfelfd(4),rfphon(4),rfstrs(4)
 real(dp) :: qphnrm(3),qphon(3,3)

! *********************************************************************

 ! Extract the block with the total energy
 qphon(:,:) = zero
 qphnrm(:) = zero
 rfphon(:) = 0
 rfelfd(:) = 0
 rfstrs(:) = 0
 rftyp = BLKTYP_d0E_xx

 call ddb%get_block(iblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,rftyp)

 if (iblok /= 0) then
   etotal = ddb%val(1,1,iblok)
 else
   etotal = huge(one)
 end if

end function ddb_get_etotal
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_get_dielt_zeff
!!
!! NAME
!!  ddb_get_dielt_zeff
!!
!! FUNCTION
!! Reads the Dielectric Tensor and the Effective Charges from the DDB file
!! Impose the charge neutrality on the effective charges and eventually select some parts of the effective charges
!!
!! INPUTS
!!  Crystal<type(crystal_t)>=Crystal structure parameters
!!  rftyp  = 1 if non-stationary block
!!           2 if stationary block
!!           3 if third order derivatives
!!  chneut=(0 => no ASR, 1 => equal repartition,2 => weighted repartition )
!!  selectz=selection of some parts of the effective charge tensor attached to one atom.
!!    (0=> no selection, 1=> trace only, 2=> symmetric part only)
!!
!! SIDE EFFECTS
!!  ddb<type(ddb_type)>=
!!    The block with the effective charges is modified if charge neutrality is imposed.
!!
!! OUTPUT
!!  iblok=Index of the block containing the data. 0 if block is not found.
!!  dielt(3,3) = Macroscopic dielectric tensor
!!  zeff(3,3,natom)=effective charge on each atom, versus electric field and atomic displacement
!!  [zeff_raw(3,3,natom)]=effective charge on each atom before enforcing charge-neutrality.
!!
!! NOTES
!!  dielt and zeff are initialized to one_3D and zero if the derivatives are not available in the DDB file.
!!
!! SOURCE

integer function ddb_get_dielt_zeff(ddb, crystal, rftyp, chneut, selectz, dielt, zeff, zeff_raw) result(iblok)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: rftyp,chneut,selectz
 class(ddb_type),intent(inout) :: ddb
 type(crystal_t),intent(in) :: crystal
!arrays
 real(dp),intent(out) :: dielt(3,3),zeff(3,3,crystal%natom)
 real(dp),optional,intent(out) :: zeff_raw(3,3,crystal%natom)

!Local variables -------------------------
!scalars
 integer :: ii
 character(len=500) :: msg
!arrays
 integer :: rfelfd(4),rfphon(4),rfstrs(4)
 real(dp) :: qphnrm(3),qphon(3,3), my_zeff_raw(3,3,crystal%natom)

! *********************************************************************

 ! Look for the Gamma Block in the DDB
 qphon(:,1)=zero
 qphnrm(1)=zero
 rfphon(1:2)=1
 rfelfd(1:2)=2
 rfstrs(1:2)=0

 !write(std_out,*)"ddb%mpert",ddb%mpert
 call ddb%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, rftyp)

 ! Compute effective charges and dielectric tensor only if the Gamma-blok was found in the DDB
 ! In case it was not found, iblok = 0
 zeff=zero; dielt=zero; dielt(1,1)=one; dielt(2,2)=one; dielt(3,3)=one
 my_zeff_raw = zero

 if (iblok /= 0) then
   write(msg, '(2a,(80a),4a)' ) ch10,('=',ii=1,80),ch10,ch10,&
   ' Dielectric Tensor and Effective Charges ',ch10
   call wrtout([std_out, ab_out], msg)

   ! Make the imaginary part of the Gamma block vanish
   write(msg, '(5a)'  ) ch10,&
   ' anaddb : Zero the imaginary part of the Dynamical Matrix at Gamma,',ch10,&
   '   and impose the ASR on the effective charges ',ch10
   call wrtout([std_out, ab_out], msg)

   ! Extrac Zeff before enforcing sum rule.
   call dtech9(ddb%val, dielt, iblok, ddb%mpert, ddb%natom, ddb%nblok, my_zeff_raw, unit=dev_null)

   ! Impose the charge neutrality on the effective charges and eventually select some parts of the effective charges
   call chneu9(chneut,ddb%val(:,:,iblok),ddb%mpert,ddb%natom,ddb%ntypat,selectz,Crystal%typat,Crystal%zion)

   ! Extraction of the dielectric tensor and the effective charges
   call dtech9(ddb%val, dielt, iblok, ddb%mpert, ddb%natom, ddb%nblok, zeff)
 end if ! iblok not found

 if (present(zeff_raw)) zeff_raw = my_zeff_raw

end function ddb_get_dielt_zeff
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_get_dielt
!!
!! NAME
!!  ddb_get_dielt
!!
!! FUNCTION
!! Reads the electronic dielectric tensor from the DDB file
!!
!! INPUTS
!!  ddb<type(ddb_type)>=Derivative database.
!!  rftyp  = 1 if non-stationary block
!!           2 if stationary block
!!           3 if third order derivatives
!!
!! OUTPUT
!!  dielt(3,3) = Macroscopic dielectric tensor (electronic contribution)
!!  iblok=Index of the block containing the data. 0 if block is not found.
!!
!! NOTES
!!  dielt is initialized to one_3D if the derivatives are not available in the DDB file.
!!
!! SOURCE

integer function ddb_get_dielt(ddb, rftyp, dielt) result(iblok)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: rftyp
 class(ddb_type),intent(in) :: ddb
!arrays
 real(dp),intent(out) :: dielt(3,3)

!Local variables -------------------------
!scalars
 integer :: mpert
 character(len=1000) :: msg
!arrays
 integer :: rfelfd(4),rfphon(4),rfstrs(4)
 real(dp) :: qphnrm(3),qphon(3,3)
 real(dp),allocatable :: tmpval(:,:,:,:)

! *********************************************************************

 ! Look for the Gamma Block in the DDB
 qphon(:,1)=zero
 qphnrm(1)=zero
 rfphon(1:2)=0
 rfelfd(1:2)=2
 rfstrs(1:2)=0

 call ddb%get_block(iblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,rftyp)

 ! Read the dielectric tensor only if the Gamma-block was found in the DDB
 ! In case it was not found, iblok = 0
 dielt=zero; dielt(1,1)=one; dielt(2,2)=one; dielt(3,3)=one

 if (iblok/=0) then
   !Extration of dielectric tensor
   mpert = ddb%mpert

   ABI_MALLOC(tmpval,(3,mpert,3,mpert))
   tmpval(:,:,:,:) = reshape(ddb%val(1,:,iblok), shape = (/3,mpert,3,mpert/))
   dielt=tmpval(1:3,ddb%natom+2,1:3,ddb%natom+2)

   write(msg,'(a,3es16.6,3es16.6,3es16.6)' )&
       ' Dielectric Tensor ',&
       dielt(1,1),dielt(1,2),dielt(1,3),&
       dielt(2,1),dielt(2,2),dielt(2,3),&
       dielt(3,1),dielt(3,2),dielt(3,3)

   call wrtout(std_out,msg)

   ABI_FREE(tmpval)
 end if ! iblok not found

end function ddb_get_dielt
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_get_quadrupoles
!!
!! NAME
!!  ddb_get_quadrupoles
!!
!! FUNCTION
!! Reads the Dynamic Quadrupoles or the P^(1) tensor from the DDB file
!!
!! INPUTS
!!  ddb<type(ddb_type)>=Derivative database.
!!  ddb_version = 6 digit integer giving date. To mantain compatibility with old DDB files.
!!  lwsym  = 0 do not symmetrize the tensor wrt efield and qvec derivative
!!             |-> 1st gradient of polarization response to atomic displacement
!!         = 1 symmetrize the tensor wrt efield and qvec derivative
!!             |-> dynamic quadrupoles
!!  rftyp  = 1 if non-stationary block
!!           2 if stationary block
!!           3 if third order derivatives
!!          33 if long wave third order derivatives
!!
!! OUTPUT
!!  quadrupoles(3,3,3,natom) = Dynamic Quadrupole tensor
!!  iblok=Index of the block containing the data. 0 if block is not found.
!!
!! NOTES
!!  quadrupoles is initialized to zero if the derivatives are not available in the DDB file.
!!
!! SOURCE

integer function ddb_get_quadrupoles(ddb, ddb_version, lwsym, rftyp, quadrupoles) result(iblok)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: ddb_version,lwsym,rftyp
 class(ddb_type),intent(in) :: ddb
!arrays
 real(dp),intent(out) :: quadrupoles(3,3,3,ddb%natom)

!Local variables -------------------------
!scalars
 character(len=500) :: msg
!arrays
 integer :: rfelfd(4),rfphon(4),rfstrs(4)
 integer :: rfqvec(4)
 real(dp) :: qphnrm(3),qphon(3,3)

! *********************************************************************

 ! Look for the Gamma Block in the DDB
 qphon(:,:)=zero
 qphnrm(:)=one
 rfphon(:)=0
 rfelfd(:)=0
 rfqvec(:)=0
 rfstrs(:)=0
 rfelfd(1)=2
 rfphon(2)=1
 rfqvec(3)=1

 call ddb%get_block(iblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,rftyp,rfqvec=rfqvec)

 ! Compute the quadrupole tensor only if the Gamma-block was found in the DDB
 ! In case it was not found, iblok = 0
 quadrupoles=zero

 if (iblok /= 0) then
   write(msg,'(2a)') ch10, ' Extract quadrupoles or P^(1) coefficients from 3DTE'
   call wrtout(std_out, msg)

   if (lwsym==1) then
     write(msg, '(3a)' ) ch10, ' Dynamical Quadrupoles Tensor (units: e Bohr)',ch10
   else if (lwsym==0) then
     write(msg, '(3a)' ) ch10, &
     ' First moment of Polarization induced by atomic displacement (1/ucvol factor not included) (units: e Bohr) ',ch10
   endif
   call wrtout([std_out, ab_out], msg)

   call dtqdrp(ddb%val(:,:,iblok),ddb_version,lwsym,ddb%mpert,ddb%natom,quadrupoles)
 end if

end function ddb_get_quadrupoles
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_get_dchidet
!!
!! NAME
!!  ddb_get_dchidet
!!
!! FUNCTION
!! Reads the non-linear optical susceptibility tensor and the
!! first-order change in the linear dielectric susceptibility
!!
!! INPUTS
!! ddb<type(ddb_type)>=Derivative Database.
!! ramansr= if /= 0, impose sum rule on first-order derivatives
!!                   of the electronic susceptibility with respect
!!                   to atomic displacements
!! nlflag= if =3, only the non-linear optical susceptibilities is computed
!!
!! OUTPUT
!! dchide(3,3,3) = non-linear optical coefficients
!! dchidt(natom,3,3,3) = first-order change of the electronic dielectric
!!   tensor induced by an individual atomic displacement
!! iblok=Index of the block containing the data. 0 if block is not found.
!!   The caller should check the returned value.
!!
!! SOURCE

integer function ddb_get_dchidet(ddb, ramansr, nlflag, dchide, dchidt) result(iblok)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: ramansr, nlflag
 class(ddb_type),intent(in) :: ddb
!arrays
 real(dp),intent(out) :: dchide(3,3,3),dchidt(ddb%natom,3,3,3)

!Local variables -------------------------
!scalars
 integer :: rftyp
!arrays
 integer :: rfelfd(4),rfphon(4),rfstrs(4)
 real(dp) :: qphnrm(3),qphon(3,3)

! *********************************************************************

 qphon(:,:) = zero
 qphnrm(:)  = one
! rfphon(1)  = 1 ; rfphon(2:3) = 0
 rfelfd(:)  = 2
 rfstrs(:)  = 0
 rftyp = BLKTYP_d3E_xx

 if (nlflag < 3) then
   rfphon(1)  = 1 ; rfphon(2:3) = 0
 else
   rfphon(1)  = 0 ; rfphon(2:3) = 0
 end if

 call ddb%get_block(iblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,rftyp)

 if (iblok /= 0) then
   call dtchi(ddb%val(:,:,iblok),dchide,dchidt,ddb%mpert,ddb%natom,ramansr,nlflag)
 else
   ! Let the caller handle the error.
   dchide = huge(one); dchidt = huge(one)
 end if

end function ddb_get_dchidet
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_get_pel
!!
!! NAME
!!  ddb_get_pel
!!
!! FUNCTION
!! Get the electronic polarizability vector from the database.
!!
!! INPUTS
!!  ddb<type(ddb_type)>=Derivative database.
!!  relaxat
!!    0 => without relaxation of the atoms 
!!    1 => with relaxation of the atoms
!!  relaxstr
!!    0 => without relaxed lattice constants
!!    1 => with relaxed lattice constants at constrained polarization
!!
!! OUTPUT
!!  pel(3) = Macroscopic polarizability vector (electronic contribution)
!!  iblok=Index of the block containing the data. 0 if block is not found.
!!
!! NOTES
!!
!! SOURCE

integer function ddb_get_pel(ddb, pel, relaxat, relaxstr) result(iblok)

!Arguments -------------------------------
!scalars
 class(ddb_type),intent(in) :: ddb
 integer, intent(in) :: relaxat
 integer, intent(in) :: relaxstr
!arrays
 real(dp),intent(out) :: pel(3)

!Local variables -------------------------
!scalars
 integer :: natom
!arrays
 integer :: rfelfd(4),rfphon(4),rfstrs(4)
 real(dp) :: qphnrm(3),qphon(3,3)

! *********************************************************************

 natom = ddb%natom

 qphon(:,:) = zero; qphnrm(:) = zero
 rfphon(:) = 0; rfstrs(:) = 0; rfelfd(:) = 2
 if (relaxat == 1) rfphon(:) = 1
 if (relaxstr == 1) rfstrs(:) = 3

 call ddb%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, BLKTYP_d1E_xx)

 if (iblok/=0) then
   pel(1:3) = ddb%val(1, 3*natom+4:3*natom+6, iblok)
 end if

end function ddb_get_pel
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_get_gred
!!
!! NAME
!!  ddb_get_gred
!!
!! FUNCTION
!! Get the forces in reduced coordinates (Hartree).
!!
!! INPUTS
!!  ddb<type(ddb_type)>=Derivative database.
!!  relaxat
!!    0 => without relaxation of the atoms 
!!    1 => with relaxation of the atoms
!!  relaxstr
!!    0 => without relaxed lattice constants
!!    1 => with relaxed lattice constants at constrained polarization
!!
!! OUTPUT
!!  gred(3,natom)=the gradient of the total energy with respect
!!                to change of reduced coordinates
!!  iblok=Index of the block containing the data. 0 if block is not found.
!!
!! NOTES
!!
!! SOURCE

integer function ddb_get_gred(ddb, gred, relaxat, relaxstr) result(iblok)

!Arguments -------------------------------
!scalars
 class(ddb_type),intent(in) :: ddb
 integer, intent(in) :: relaxat
 integer, intent(in) :: relaxstr
!arrays
 real(dp),intent(out),allocatable :: gred(:,:)

!Local variables -------------------------
!scalars
 integer :: natom
 integer :: idir, iatom, index
!arrays
 integer :: rfelfd(4),rfphon(4),rfstrs(4)
 real(dp) :: qphnrm(3),qphon(3,3)

! *********************************************************************

 natom = ddb%natom

 ABI_MALLOC(gred, (3, natom))

 qphon(:,:) = zero; qphnrm(:) = zero
 rfphon(:) = 0; rfstrs(:) = 0; rfelfd(:) = 2
 if (relaxat == 1) rfphon(:) = 1
 if (relaxstr == 1) rfstrs(:) = 3

 ! GA: I dont see why relaxstr is relevant as an input
 call ddb%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, BLKTYP_d1E_xx)

 if (iblok/=0) then
   if (relaxat == 1) then
     index = 0
     do iatom = 1, natom
       do idir = 1, 3
         index = index+1
         gred(idir, iatom) = ddb%val(1, index, iblok)
       end do
     end do
   end if
 end if

end function ddb_get_gred
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_get_strten
!!
!! NAME
!!  ddb_get_strten
!!
!! FUNCTION
!! Get the stress tensor.
!!
!! INPUTS
!!  ddb<type(ddb_type)>=Derivative database.
!!  relaxat
!!    0 => without relaxation of the atoms 
!!    1 => with relaxation of the atoms
!!  relaxstr
!!    0 => without relaxed lattice constants
!!    1 => with relaxed lattice constants at constrained polarization
!!
!! OUTPUT
!!  strten(6)=the stress tensor in cartesian coordinates.
!!  iblok=Index of the block containing the data. 0 if block is not found.
!!
!! NOTES
!!
!! SOURCE

integer function ddb_get_strten(ddb, strten, relaxat, relaxstr) result(iblok)

!Arguments -------------------------------
!scalars
 class(ddb_type),intent(in) :: ddb
 integer, intent(in) :: relaxat
 integer, intent(in) :: relaxstr
!arrays
 real(dp),intent(out) :: strten(6)

!Local variables -------------------------
!scalars
 integer :: natom
 integer :: ii, index
!arrays
 integer :: rfelfd(4),rfphon(4),rfstrs(4)
 real(dp) :: qphnrm(3),qphon(3,3)

! *********************************************************************

 natom = ddb%natom

 qphon(:,:) = zero; qphnrm(:) = zero
 rfphon(:) = 0; rfstrs(:) = 0; rfelfd(:) = 2
 if (relaxat == 1) rfphon(:) = 1
 if (relaxstr == 1) rfstrs(:) = 3

 ! GA: I dont see why relaxstr is relevant as an input
 call ddb%get_block(iblok, qphon, qphnrm, rfphon, rfelfd, rfstrs, BLKTYP_d1E_xx)

 if (iblok/=0) then
   if (relaxstr == 1) then
     index = 3*natom+6
     do ii = 1, 6
       index = index+1
       strten(ii) = ddb%val(1, index, iblok)
     end do
   end if
 end if

end function ddb_get_strten
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_get_asrq0
!! NAME
!!  ddb_get_asrq0
!!
!! FUNCTION
!!  In case the interatomic forces are not calculated, the
!!  ASR-correction has to be determined here from the Dynamical matrix at Gamma.
!!  In case the DDB does not contain this information, the subroutine returns iblok=0
!!  %d2asr is initialized and set to zero to preserve the old behaviour.
!!
!! INPUTS
!!  asr=Input variable selecting the method for the ASR
!!  rftyp  = 1 if non-stationary block
!!           2 if stationary block
!!           3 if third order derivatives
!!  xcart(3,ddb%atom)=Cartesian coordinates of the atoms.
!!
!! SIDE EFFECTS
!!  ddb<type(ddb_type)>= Database with the derivates. The routine does not change it
!!  except when asr is in [3,4]. TODO This should not happen.
!!
!! OUTPUT
!! asrq0<asrq0_t>
!!   iblok= is set to 0 if the Gamma block is not found
!!
!! SOURCE

type(asrq0_t) function ddb_get_asrq0(ddb, asr, rftyp, xcart) result(asrq0)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: asr,rftyp
 class(ddb_type),intent(inout) :: ddb
!arrays
 real(dp),intent(in) :: xcart(3,ddb%natom)

!Local variables-------------------------------
!scalars
 integer :: dims,iblok
 !character(len=500) :: msg
!arrays
 integer :: rfelfd(4),rfphon(4),rfstrs(4)
 real(dp) :: qphnrm(3),qphon(3,3)
 real(dp),allocatable :: d2asr_res(:,:,:,:,:),d2cart(:,:)

! ************************************************************************

 asrq0%asr = asr; asrq0%natom = ddb%natom

 ! Find the Gamma block in the DDB (no need for E-field entries)
 qphon(:,1)=zero
 qphnrm(1)=zero
 rfphon(1:2)=1
 rfelfd(:)=0
 rfstrs(:)=0

 call ddb%get_block(asrq0%iblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,rftyp)
 ! this is to maintain the old behaviour in which the arrays where allocated and set to zero in anaddb.
 ABI_MALLOC(asrq0%d2asr, (2,3,ddb%natom,3,ddb%natom))
 asrq0%d2asr = zero

 if (asrq0%iblok == 0) return
 iblok = asrq0%iblok

 select case (asrq0%asr)
 case (0)
   continue

 case (1,2)
   call asria_calc(asr,asrq0%d2asr,ddb%val(:,:,iblok),ddb%mpert,ddb%natom)

 case (3,4)
   ! Rotational invariance for 1D and 0D systems
   ! Compute uinvers, vtinvers and singular matrices.
   dims = 3*ddb%natom*(3*ddb%natom-1) / 2
   ABI_CALLOC(asrq0%uinvers, (dims, dims))
   ABI_CALLOC(asrq0%vtinvers,(dims, dims))
   ABI_CALLOC(asrq0%singular, (dims))

   call asrprs(asr,1,3,asrq0%uinvers,asrq0%vtinvers,asrq0%singular,&
     ddb%val(:,:,iblok),ddb%mpert,ddb%natom,xcart)

 case (5)
   ! d2cart is a temp variable here
   ABI_MALLOC(d2cart,(2,ddb%msize))
   d2cart = ddb%val(:,:,iblok)
   ! calculate diagonal correction
   call asria_calc(2,asrq0%d2asr,d2cart,ddb%mpert,ddb%natom)
   ! apply diagonal correction
   call asria_corr(2,asrq0%d2asr,d2cart,ddb%mpert,ddb%natom)
   ! hermitianize
   call mkherm(d2cart,3*ddb%mpert)
   ! remove remaining ASR rupture due to Hermitianization
   ABI_MALLOC(d2asr_res,(2,3,ddb%natom,3,ddb%natom))
   call asria_calc(asr,d2asr_res,d2cart,ddb%mpert,ddb%natom)
   ! full correction is sum of both
   asrq0%d2asr = asrq0%d2asr + d2asr_res

   ABI_FREE(d2cart)
   ABI_FREE(d2asr_res)

 case default
   ABI_ERROR(sjoin("Wrong value for asr:", itoa(asr)))
 end select

end function ddb_get_asrq0
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_symmetrize_and_transform
!! NAME
!!  ddb_symmetrize_and_transform
!!
!! FUNCTION
!! First apply symmetry operations, then
!! transform the second-derivative matrix from reduced
!! coordinates to cartesian coordinates, and also
!! 1) add the ionic part of the effective charges,
!! 2) normalize the electronic dielectric tensor, and
!!    add the vacuum polarisation
!!
!! INPUTS
!!  crystal<type(crystal_t)>
!!  iblock=the block index on which to act
!!
!! SIDE EFFECTS
!!  ddb<type(ddb_type)>= 
!!
!! OUTPUT
!!
!! SOURCE

subroutine ddb_symmetrize_and_transform(ddb, crystal, iblok)

!Arguments -------------------------------
!scalars
 class(ddb_type),intent(inout) :: ddb
 class(crystal_t),intent(in) :: crystal
 integer,intent(in) :: iblok
 !integer,intent(inout) :: indsym(4,msym,natom)
 !integer,intent(out) :: symrec(3,3,msym),symrel(3,3,msym),symafm(msym)

!Local variables-------------------------------
!scalars
 integer :: mpert,natom,nsym,ntypat
 integer :: nsize,timrev
 integer :: i1dir,i1pert,i2dir,i2pert,i3dir,i3pert
!arrays
 !integer :: symq(4,2,msym)
 integer,allocatable :: symq(:,:,:)
 integer,allocatable :: car3flg(:,:,:,:,:,:),carflg(:,:,:,:)
 integer,allocatable :: tmpflg(:,:,:,:,:,:),rfpert(:,:,:,:,:,:)
 real(dp) :: gprimd(3,3),qpt(3),rprimd(3,3)
 real(dp),allocatable :: d2cart(:,:,:,:,:),d3cart(:,:,:,:,:,:,:)
 real(dp),allocatable :: tmpval(:,:,:,:,:,:,:)

! ************************************************************************

 mpert = ddb%mpert
 natom = ddb%natom
 ntypat = crystal%ntypat

 nsym = crystal%nsym
 rprimd(:,:) = crystal%rprimd(:,:)
 gprimd(:,:) = crystal%gprimd(:,:)

 ABI_MALLOC(symq,(4,2,nsym))

 !  Here complete the matrix by symmetrisation of the existing elements
 if (is_type_d2E(ddb%typ(iblok))) then

   qpt(1)=ddb%qpt(1,iblok)/ddb%nrm(1,iblok)
   qpt(2)=ddb%qpt(2,iblok)/ddb%nrm(1,iblok)
   qpt(3)=ddb%qpt(3,iblok)/ddb%nrm(1,iblok)

   ! Examine the symmetries of the q wavevector
   call littlegroup_q(crystal%nsym,qpt,symq,crystal%symrec,crystal%symafm,timrev,prtvol=0)

   !GA: Note that d2sym3 and cart29 expect different shapes for tmpflg and tmpval
   !    hence the extra dimensions
   nsize=3*mpert*3*mpert
   ABI_MALLOC(tmpflg,(3,mpert,3,mpert,1,1))
   ABI_MALLOC(tmpval,(2,3,mpert,3,mpert,1,1))

   tmpflg(:,:,:,:,1,1) = reshape(ddb%flg(1:nsize,iblok), shape = (/3,mpert,3,mpert/))
   tmpval(1,:,:,:,:,1,1) = reshape(ddb%val(1,1:nsize,iblok), shape = (/3,mpert,3,mpert/))
   tmpval(2,:,:,:,:,1,1) = reshape(ddb%val(2,1:nsize,iblok), shape = (/3,mpert,3,mpert/))

   ! Then apply symmetry operations
   call d2sym3(tmpflg,tmpval,crystal%indsym,mpert,natom,nsym,qpt,symq,crystal%symrec,crystal%symrel,timrev,1)

   ! Transform the dynamical matrix in cartesian coordinates
   ABI_MALLOC(carflg,(3,mpert,3,mpert))
   ABI_MALLOC(d2cart,(2,3,mpert,3,mpert))

   call cart29(tmpflg,tmpval,carflg,d2cart,gprimd,1,mpert,natom,1,ntypat,rprimd,crystal%typat,crystal%ucvol,crystal%zion)

   ddb%flg(1:nsize,iblok) = reshape(carflg,shape = (/3*mpert*3*mpert/))
   ddb%val(1,1:nsize,iblok) = reshape(d2cart(1,:,:,:,:), shape = (/3*mpert*3*mpert/))
   ddb%val(2,1:nsize,iblok) = reshape(d2cart(2,:,:,:,:), shape = (/3*mpert*3*mpert/))

   ABI_FREE(carflg)
   ABI_FREE(d2cart)
   ABI_FREE(tmpflg)
   ABI_FREE(tmpval)

 else if (ddb%typ(iblok) == BLKTYP_d3E_xx) then

   nsize=3*mpert*3*mpert*3*mpert
   ABI_MALLOC(tmpflg,(3,mpert,3,mpert,3,mpert))
   ABI_MALLOC(tmpval,(2,3,mpert,3,mpert,3,mpert))
   ABI_MALLOC(rfpert,(3,mpert,3,mpert,3,mpert))

   tmpflg(:,:,:,:,:,:) = reshape(ddb%flg(1:nsize,iblok), shape = (/3,mpert,3,mpert,3,mpert/))
   tmpval(1,:,:,:,:,:,:) = reshape(ddb%val(1,1:nsize,iblok), shape = (/3,mpert,3,mpert,3,mpert/))
   tmpval(2,:,:,:,:,:,:) = reshape(ddb%val(2,1:nsize,iblok), shape = (/3,mpert,3,mpert,3,mpert/))

   ! Set the elements that are zero by symmetry for raman and
   ! non-linear optical susceptibility tensors
   rfpert = 0
   rfpert(:,natom+2,:,natom+2,:,natom+2) = 1
   rfpert(:,1:natom,:,natom+2,:,natom+2) = 1
   rfpert(:,natom+2,:,1:natom,:,natom+2) = 1
   rfpert(:,natom+2,:,natom+2,:,1:natom) = 1
   call sytens(crystal%indsym,mpert,natom,nsym,rfpert,crystal%symrec,crystal%symrel)
   do i1pert = 1,mpert
     do i2pert = 1,mpert
       do i3pert = 1,mpert
         do i1dir=1,3
           do i2dir=1,3
             do i3dir=1,3
               if ((rfpert(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)==-2) .and. &
                   (tmpflg(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)/=1)) then
                 tmpval(:,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = zero
                 tmpflg(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)=1
               end if
             end do
           end do
         end do
       end do
     end do
   end do

   call d3sym(tmpflg,tmpval,crystal%indsym,mpert,natom,nsym,crystal%symrec,crystal%symrel)

   ABI_MALLOC(d3cart,(2,3,mpert,3,mpert,3,mpert))
   ABI_MALLOC(car3flg,(3,mpert,3,mpert,3,mpert))

   call nlopt(tmpflg,car3flg,tmpval,d3cart,gprimd,mpert,natom,rprimd,crystal%ucvol)

   ddb%flg(1:nsize,iblok) = reshape(car3flg, shape = (/3*mpert*3*mpert*3*mpert/))
   ddb%val(1,1:nsize,iblok) = reshape(d3cart(1,:,:,:,:,:,:), shape = (/3*mpert*3*mpert*3*mpert/))
   ddb%val(2,1:nsize,iblok) = reshape(d3cart(2,:,:,:,:,:,:), shape = (/3*mpert*3*mpert*3*mpert/))

   ABI_FREE(d3cart)
   ABI_FREE(car3flg)
   ABI_FREE(tmpflg)
   ABI_FREE(tmpval)
   ABI_FREE(rfpert)

 else if (ddb%typ(iblok) == BLKTYP_d3E_lw) then

   nsize=3*mpert*3*mpert*3*mpert
   ABI_MALLOC(tmpflg,(3,mpert,3,mpert,3,mpert))
   ABI_MALLOC(tmpval,(2,3,mpert,3,mpert,3,mpert))

   tmpflg(:,:,:,:,:,:) = reshape(ddb%flg(1:nsize,iblok), shape = (/3,mpert,3,mpert,3,mpert/))
   tmpval(1,:,:,:,:,:,:) = reshape(ddb%val(1,1:nsize,iblok), shape = (/3,mpert,3,mpert,3,mpert/))
   tmpval(2,:,:,:,:,:,:) = reshape(ddb%val(2,1:nsize,iblok), shape = (/3,mpert,3,mpert,3,mpert/))

   ABI_MALLOC(d3cart,(2,3,mpert,3,mpert,3,mpert))
   ABI_MALLOC(car3flg,(3,mpert,3,mpert,3,mpert))

   call lwcart(tmpflg,car3flg,tmpval,d3cart,gprimd,mpert,natom,rprimd)

   ddb%flg(1:nsize,iblok) = reshape(car3flg, shape = (/3*mpert*3*mpert*3*mpert/))
   ddb%val(1,1:nsize,iblok) = reshape(d3cart(1,:,:,:,:,:,:), shape = (/3*mpert*3*mpert*3*mpert/))
   ddb%val(2,1:nsize,iblok) = reshape(d3cart(2,:,:,:,:,:,:), shape = (/3*mpert*3*mpert*3*mpert/))

   ABI_FREE(d3cart)
   ABI_FREE(car3flg)
   ABI_FREE(tmpflg)
   ABI_FREE(tmpval)
 end if

 ABI_FREE(symq)

end subroutine ddb_symmetrize_and_transform
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_diagoq
!! NAME
!!  ddb_diagoq
!!
!! FUNCTION
!!  Compute the phonon frequencies at the specified q-point by performing
!!  a direct diagonalization of the dynamical matrix. The q-point **MUST** be
!!  one the points stored in the DDB file.
!!
!! INPUTS
!!  ddb<type(ddb_type)>=Object storing the DDB results.
!!  crystal<type(crystal_t)> = Information on the crystalline structure.
!!  asrq0<asrq0_t>=Object for the treatment of the ASR based on the q=0 block found in the DDB file.
!!  symdynmat=If equal to 1, the dynamical matrix is symmetrized in dfpt_phfrq before the diagonalization.
!!  rftyp  = 1 if non-stationary block
!!           2 if stationary block
!!           3 if third order derivatives
!!  qpt(3)=q-point in reduced coordinates.
!!
!! OUTPUT
!!  phfrq(3*crystal%natom)=Phonon frequencies in Hartree
!!  displ_cart(2,3*%natom,3*%natom)=Phonon displacement in Cartesian coordinates
!!  [out_eigvec(2*3*natom*3*natom) = The igenvectors of the dynamical matrix.
!!  [out_displ_red(2*3*natom*3*natom) = The displacement in reduced coordinates.
!!
!! SOURCE

subroutine ddb_diagoq(ddb, crystal, qpt, asrq0, symdynmat, rftyp, phfrq, displ_cart, &
                      out_eigvec,out_displ_red)   ! Optional [out]

!Arguments -------------------------------
!scalars
 integer,intent(in) :: symdynmat
 integer,intent(in) :: rftyp
 class(ddb_type),intent(in) :: ddb
 type(asrq0_t),intent(inout) :: asrq0
 type(crystal_t),intent(in) :: crystal
!arrays
 real(dp),intent(in) :: qpt(3)
 real(dp),intent(out) :: displ_cart(2,3,crystal%natom,3,crystal%natom)
 real(dp),intent(out) :: phfrq(3*crystal%natom)
 !real(dp),optional,intent(out) :: out_d2cart(2,3*crystal%natom,3*crystal%natom)
 real(dp),optional,intent(out) :: out_eigvec(2,3,crystal%natom,3*crystal%natom)
 real(dp),optional,intent(out) :: out_displ_red(2,3,crystal%natom,3*crystal%natom)

!Local variables-------------------------------
 integer :: iblok,natom
!arrays
 integer :: rfphon(4),rfelfd(4),rfstrs(4)
 real(dp) :: qphnrm(3), qphon_padded(3,3),d2cart(2,ddb%msize),my_qpt(3)
 real(dp) :: eigvec(2,3,crystal%natom,3*crystal%natom),eigval(3*crystal%natom)

! ************************************************************************

 ! Use my_qpt because dfpt_phfrq can change the q-point (very bad design)
 qphnrm = one; my_qpt = qpt

 ! Look for the information in the DDB (no interpolation here!)
 rfphon(1:2)=1; rfelfd(1:2)=0; rfstrs(1:2)=0
 qphon_padded = zero; qphon_padded(:,1) = qpt
 natom = crystal%natom

 call ddb%get_block(iblok, qphon_padded, qphnrm, rfphon, rfelfd, rfstrs, rftyp)
 ABI_CHECK(iblok /= 0, sjoin("Cannot find q-point ", ktoa(qpt)," in DDB file"))

 ! Copy the dynamical matrix in d2cart
 d2cart(:,1:ddb%msize) = ddb%val(:,:,iblok)

 ! Eventually impose the acoustic sum rule based on previously calculated d2asr
 call asrq0%apply(natom, ddb%mpert, ddb%msize, crystal%xcart, d2cart)

 ! Calculation of the eigenvectors and eigenvalues of the dynamical matrix
 call dfpt_phfrq(ddb%amu,displ_cart,d2cart,eigval,eigvec,crystal%indsym,&
   ddb%mpert,crystal%nsym,natom,crystal%nsym,crystal%ntypat,phfrq,qphnrm(1),my_qpt,&
   crystal%rprimd,symdynmat,crystal%symrel,crystal%symafm,crystal%typat,crystal%ucvol)

 ! Return the dynamical matrix and the eigenvector for this q-point
 !if (present(out_d2cart)) out_d2cart = d2cart(:,:3*natom,:3*natom)
 if (present(out_eigvec)) out_eigvec = eigvec

 ! Return phonon displacement in reduced coordinates.
 if (present(out_displ_red)) call phdispl_cart2red(natom, crystal%gprimd, displ_cart, out_displ_red)

end subroutine ddb_diagoq
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/asrq0_apply
!! NAME
!! asrq0_apply
!!
!! FUNCTION
!!  Impose the acoustic sum rule based on the q=0 block found in the DDB file.
!!
!! INPUTS
!!  asrq0<asrq0_t>=Object for the treatment of the ASR based on the q=0 block found in the DDB file.
!!  natom=Number of atoms per unit cell.
!!  mpert=Maximum number of perturbation (reported in ddb%mpert)
!!  msize=Maximum size of array ddb%val
!!  xcart(3,natom)=Atomic positions in Cartesian coordinates
!!
!! SIDE EFFECTS
!!   d2cart=matrix of second derivatives of total energy, in cartesian coordinates
!!   Input: Values stored in ddb%
!!   Output: Changed to enforce ASR.
!!
!! SOURCE

subroutine asrq0_apply(asrq0, natom, mpert, msize, xcart, d2cart)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: natom, msize, mpert
 class(asrq0_t),intent(inout) :: asrq0
!arrays
 real(dp),intent(in) :: xcart(3,natom)
 real(dp),intent(inout) :: d2cart(2,msize)

! ************************************************************************

 if (asrq0%asr /= 0 .and. asrq0%iblok == 0) then
   ABI_WARNING("asr != 0 but DDB file does not contain q=Gamma. D(q) cannot be corrected")
   return
 end if

 select case (asrq0%asr)
 case (0)
   return
 case (1,2,5)
   call asria_corr(asrq0%asr, asrq0%d2asr, d2cart, mpert, natom)
 case (3,4)
   ! Impose acoustic sum rule plus rotational symmetry for 0D and 1D systems
   call asrprs(asrq0%asr,2,3,asrq0%uinvers,asrq0%vtinvers,asrq0%singular,d2cart,mpert,natom,xcart)
 case default
   ABI_ERROR(sjoin("Wrong value for asr:", itoa(asrq0%asr)))
 end select

end subroutine asrq0_apply
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/asrq0_free
!! NAME
!! asrq0_free
!!
!! FUNCTION
!!   Free dynamic memory
!!
!! SOURCE

subroutine asrq0_free(asrq0)

!Arguments -------------------------------
 class(asrq0_t),intent(inout) :: asrq0

! ************************************************************************

 ! real
 ABI_SFREE(asrq0%d2asr)
 ABI_SFREE(asrq0%singular)
 ABI_SFREE(asrq0%uinvers)
 ABI_SFREE(asrq0%vtinvers)

end subroutine asrq0_free
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_write_block_txt
!!
!! NAME
!! ddb_write_block_txt
!!
!! FUNCTION
!! This routine writes blocks of data in the DDB in text format.
!!
!! INPUTS
!! choice= (2 => write), (3 => write minimal info )
!! mpert =maximum number of ipert
!! msize=maximum size of the arrays flags and values
!! nunit=unit number for the data block file
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output
!! ddb = ddb block datastructure
!! ddb%typ=type of the block:
!!   0 => total energy
!!   1 => second-order energy derivatives, non-stationary block
!!   2 => second-order energy derivatives, stationary block
!!   3 => third-order energy derivatives
!!   4 => first-order energy derivatives: forces, stresses and polarization
!!   5 => second-order eigenvalue derivatives
!! ddb%flg(msize)=flag for every matrix element (0=> the element is
!!  not in the data block), (1=> the element is in the data blok)
!! ddb%qpt(9)=wavevector of the perturbation(s). The elements from
!!  1 to 3 are used if we are dealing with the 2nd derivative of
!!  total energy (only one wavevector), while all elements are
!!  used in case of a third order derivative of total energy
!!  (three wavevector could be present)
!! ddb%nrm(3)=normalization factors for the three allowed wavevectors.
!! ddb%val(2,msize)=real(dp), complex, value of the
!!  matrix elements that are present in the data block
!! blkval2(2,msize,mband,nkpt) = value of the matrix elements
!!  that are present in a block of EIGR2D/EIGI2D
!!
!! NOTES
!! only executed by one processor.
!!
!! SOURCE

subroutine ddb_write_block_txt(ddb,iblok,choice,mband,mpert,msize,nkpt,nunit,&
                           blkval2,kpt) !optional

!Arguments -------------------------------
!scalars
 integer,intent(in) :: choice,mband,mpert,msize,nkpt,nunit
 integer,intent(in) :: iblok
 class(ddb_type),intent(in) :: ddb
!arrays
 real(dp),intent(in),optional :: kpt(3,nkpt)
 real(dp),intent(in),optional :: blkval2(2,msize,mband,nkpt)

!Local variables -------------------------
!scalars
 integer :: iband,idir1,idir2,idir3,ii,ikpt,ipert1,ipert2,ipert3
 integer :: nelmts
 logical :: eig2d_

! *********************************************************************

 ! GA: Remove choice option.
 ! choice=3 is used to write summary info at the end of the DDB.
 ! With MG and MV, we agreed that it could be removed.

 ! GA: Remove arguments: mband, blkval2, kpt
 ! This feature of writing eigenvalues 2nd deriv is no longer used
 ! (see 80_tdep/m_tdep_abitypes.F90)
 ! With FB, we agreed that it could be removed.

 eig2d_ = .false.
 if(present(blkval2).and.present(kpt)) eig2d_ = .true.


 ! Count the number of elements
 nelmts=0
 do ii=1,msize
   if(ddb%flg(ii,iblok)==1)nelmts=nelmts+1
 end do

 ! Write the block type and number of elements
 write(nunit,*)' '
 if (ddb%typ(iblok) == BLKTYP_d0E_xx) then
   write(nunit, '(a,i8)' )' Total energy                 - # elements :',nelmts
 else if (ddb%typ(iblok)==BLKTYP_d2E_ns) then
   write(nunit, '(a,i8)' )' 2nd derivatives (non-stat.)  - # elements :',nelmts
 else if(ddb%typ(iblok)==BLKTYP_d2E_st) then
   write(nunit, '(a,i8)' )' 2nd derivatives (stationary) - # elements :',nelmts
 else if(ddb%typ(iblok)==BLKTYP_d3E_xx) then
   write(nunit, '(a,i8)' )' 3rd derivatives              - # elements :',nelmts
 else if (ddb%typ(iblok) == BLKTYP_d1E_xx) then
   write(nunit, '(a,i8)' )' 1st derivatives              - # elements :',nelmts
 else if (ddb%typ(iblok) == BLKTYP_d2eig_re) then
   write(nunit, '(a,i8)' )' 2nd eigenvalue derivatives   - # elements :',nelmts
 else if(ddb%typ(iblok)==BLKTYP_d3E_lw) then
   write(nunit, '(a,i8)' )' 3rd derivatives (long wave)  - # elements :',nelmts
 end if

 ! Write the 2nd derivative block
 if (is_type_d2E(ddb%typ(iblok))) then

   ! Write the phonon wavevector
   write(nunit, '(a,3es16.8,f6.1)' )' qpt',(ddb%qpt(ii,iblok),ii=1,3),ddb%nrm(1,iblok)

   ! Write the matrix elements
   if(choice==2)then
     ii=0
     do ipert2=1,mpert
       do idir2=1,3
         do ipert1=1,mpert
           do idir1=1,3
             ii=ii+1
             if(ddb%flg(ii,iblok)==1)then
               write(nunit,'(4i4,2d22.14)') idir1, ipert1, idir2, ipert2, &
                ddb%val(1,ii,iblok),ddb%val(2,ii,iblok)
             end if
           end do
         end do
       end do
     end do
   end if


 else if (is_type_d3E(ddb%typ(iblok))) then
   ! Write the 3rd derivative block

   ! Write the phonon wavevectors
   write(nunit, '(a,3es16.8,f6.1)' )' qpt',(ddb%qpt(ii,iblok),ii=1,3),ddb%nrm(1,iblok)
   write(nunit, '(a,3es16.8,f6.1)' )'    ',(ddb%qpt(ii,iblok),ii=4,6),ddb%nrm(2,iblok)
   write(nunit, '(a,3es16.8,f6.1)' )'    ',(ddb%qpt(ii,iblok),ii=7,9),ddb%nrm(3,iblok)

   ! Write the matrix elements
   if(choice==2)then
     ii=0
     do ipert3=1,mpert
       do idir3=1,3
         do ipert2=1,mpert
           do idir2=1,3
             do ipert1=1,mpert
               do idir1=1,3
                 ii=ii+1
                 if(ddb%flg(ii,iblok)==1)then
                   write(nunit, '(6i4,2d22.14)' )&
                    idir1,ipert1,idir2,ipert2,idir3,ipert3,ddb%val(1,ii,iblok),ddb%val(2,ii,iblok)
                 end if
               end do
             end do
           end do
         end do
       end do
     end do
   end if


 else if (is_type_d0E(ddb%typ(iblok))) then
   !  Write total energy
   if (choice == 2) write(nunit,'(2d22.14)')ddb%val(1,1,iblok),ddb%val(2,1,iblok)

 else if (is_type_d1E(ddb%typ(iblok))) then
   !  Write the 1st derivative blok
   if (choice == 2) then
     ii = 0
     do ipert1 = 1, mpert
       do idir1 = 1, 3
         ii = ii + 1
         if (ddb%flg(ii,iblok) == 1) then
           write(nunit,'(2i4,2d22.14)')idir1,ipert1,ddb%val(1,ii,iblok),ddb%val(2,ii,iblok)
         end if
       end do
     end do
   end if

 else if (is_type_d2eig(ddb%typ(iblok))) then
   ! Write the phonon wavevector
   write(nunit, '(a,3es16.8,f6.1)' )' qpt',(ddb%qpt(ii,iblok),ii=1,3),ddb%nrm(1,iblok)
   ! Write the matrix elements
   ! GA: Note that isppol is invisible here. It is simply marked as more bands.
   !     To be changed in a future version of text format.
   if(choice==2)then
     if (eig2d_) then
       do ikpt=1,nkpt
         write(nunit,'(a,3es16.8)')' K-point:',(kpt(ii,ikpt),ii=1,3)
         do iband=1,mband
           write(nunit,'(a,i3)')' Band:',iband
           ii=0
           do ipert2=1,mpert
             do idir2=1,3
               do ipert1=1,mpert
                 do idir1=1,3
                   ii=ii+1
                   if(ddb%flg(ii,iblok)==1)then
                     write(nunit,'(4i4,2d22.14)')idir1,ipert1,idir2,ipert2,blkval2(1,ii,iband,ikpt),blkval2(2,ii,iband,ikpt)
                   end if
                 end do !idir1
               end do !ipert1
             end do !idir2
           end do !ipert2
         end do !iband
       end do !ikpt
     end if !eig2d_
   end if !choice
 end if !ddb%typ(iblok)

end subroutine ddb_write_block_txt
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_write
!! NAME
!! ddb_write
!!
!! FUNCTION
!!  Write the DDB file in either txt or netcdf format.
!!
!! INPUTS
!!  ddb_hdr=ddb header object.
!!  filename=name of the file being written (abo_DS*_DDB)
!!  with_psps
!!      1-> include information on pseudopoentials
!!      0-> do not include information on pseudopoentials
!!  comm=MPI communicator
!!
!! SOURCE

subroutine ddb_write(ddb, ddb_hdr, filename, with_psps, comm)

!Arguments -------------------------------
 class(ddb_type),intent(inout) :: ddb
 type(ddb_hdr_type),intent(inout) :: ddb_hdr
 character(len=fnlen),intent(in) :: filename
 integer,intent(in),optional :: with_psps
 integer,intent(in),optional :: comm

!Local variables-------------------------------
 character(len=fnlen) :: filename_
 integer :: iomode

! ************************************************************************

  call ddb_hdr%get_iomode(filename, 2, iomode, filename_)

  if (iomode==IO_MODE_ETSF) then
    call ddb%write_nc(ddb_hdr, filename_, comm=comm, with_psps=with_psps)
  else if (iomode==IO_MODE_FORTRAN) then
    call ddb%write_txt(ddb_hdr, filename_, with_psps=with_psps, comm=comm)
    ddb_hdr%mpert = ddb%mpert  ! Text format doesnt know about mpert.
  end if

end subroutine ddb_write
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_write_txt
!! NAME
!! ddb_write_txt
!!
!! FUNCTION
!!  Write the DDB file in text format.
!!
!! INPUTS
!!  ddb_hdr=ddb header object.
!!  filename=name of the file being written (abo_DS*_DDB)
!!  with_psps
!!      1-> include information on pseudopoentials
!!      0-> do not include information on pseudopoentials
!!
!! SOURCE

subroutine ddb_write_txt(ddb, ddb_hdr, filename, with_psps, comm)

!Arguments -------------------------------
 class(ddb_type),intent(inout) :: ddb
 type(ddb_hdr_type),intent(inout) :: ddb_hdr
 character(len=*),intent(in) :: filename
 integer,intent(in),optional :: with_psps
 integer,intent(in),optional :: comm

!Local variables -------------------------
!scalars
 integer :: iblok
 integer,parameter :: master=0, choice=2

! ************************************************************************

  if (present(comm)) then
    if (xmpi_comm_rank(comm) /= master) return
  end if

 call ddb_hdr%open_write_txt(filename, with_psps)

 do iblok=1,ddb%nblok
   call ddb%write_block_txt(iblok,choice,1,ddb%mpert,ddb%msize,ddb_hdr%nkpt,ddb_hdr%unddb)
 end do

 call ddb_hdr%close()

end subroutine ddb_write_txt
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_write_d2eig
!! NAME
!! ddb_write_d2eig
!!
!! FUNCTION
!!  Write the current eig2d data as the next block in the ddb file.
!!
!! INPUTS
!!  ddb_hdr=ddb header object.
!!  unddb=unit of the open ddb file in text format or netcdf identifier.
!!
!!
!! SOURCE

subroutine ddb_write_d2eig(ddb, ddb_hdr, iblok, comm)
!Arguments -------------------------------
 class(ddb_type),intent(inout) :: ddb
 type(ddb_hdr_type),intent(inout) :: ddb_hdr
 integer,intent(in) :: iblok
 integer,intent(in),optional :: comm

!Local variables -------------------------
!scalars
 integer,parameter :: master=0
 character(len=500) :: msg

! ************************************************************************

  if (present(comm)) then
    if (xmpi_comm_rank(comm) /= master) return
  end if

  if (ddb_hdr%has_open_file_nc) then

    call ddb%write_d2eig_nc(ddb_hdr%ncid, iblok)

  else if (ddb_hdr%has_open_file_txt) then

    call ddb%write_d2eig_txt(ddb_hdr%unddb, iblok)

  else 
    write(msg, '(3a)' )&
    ! File has not been opened by ddb_hdr
    'Attempting to write into unopen DDB file.',ch10,&
    'Action: contact Abinit group.'
    ABI_ERROR(msg)
  end if

end subroutine ddb_write_d2eig
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_write_d2eig_nc
!! NAME
!! ddb_write_d2eig_nc
!!
!! FUNCTION
!!  Write the current d2eig data in the ddb netcdf file.
!!
!! INPUTS
!!  iblok=index of the eig2d block within the d2eig subgroup.
!!  ncid=netcdf identifier of a file open in writing mode.
!!  comm=MPI communicator.
!!
!! SOURCE

subroutine ddb_write_d2eig_nc(ddb, ncid, iblok, comm)
!Arguments -------------------------------
 class(ddb_type),intent(inout) :: ddb
 integer,intent(in) :: ncid
 integer,intent(in) :: iblok
 integer,intent(in),optional :: comm

!Local variables -------------------------
!scalars
 integer,parameter :: master=0
 integer :: iband, jband, bandshift, isppol, mband
 integer :: ikpt, ipert1, idir1, ipert2, idir2, ii
 integer :: ncid_d2eig, ncerr
 real(dp) :: qpt(3)
 real(dp), allocatable :: matrix_d2eig(:,:,:,:,:,:,:)
 real(dp), allocatable :: matrix_d2eig_isppol(:,:,:,:,:,:,:)
 integer, allocatable :: flg_d2eig(:,:,:,:)

! ************************************************************************


  if (present(comm)) then
    if (xmpi_comm_rank(comm) /= master) return
  end if

  ncid_d2eig = nctk_idgroup(ncid, 'd2eig')

  qpt(1:3) = ddb%qpt(1:3,iblok)
  ncerr = nf90_put_var(ncid_d2eig, nctk_idname(ncid_d2eig,&
                         'reduced_coordinates_of_qpoints'),&
                         qpt,&
                         start=[1,iblok])
  NCF_CHECK(ncerr)
  ncerr = nf90_put_var(ncid_d2eig, nctk_idname(ncid_d2eig,&
                         'qpoints_normalization'),&
                         ddb%nrm(1,iblok),&
                         start=[iblok])
  NCF_CHECK(ncerr)

  ! GA: Here we assume that all blocks are d2eig blocks.
  !     Otherwise, iblok is only the 'local' iblok index
  !     and we would need to figure out the corresponding 'global' index
  call ddb%get_d2eig(matrix_d2eig, flg_d2eig, iblok)

  mband = ddb%nband / ddb%nsppol
  ABI_MALLOC(matrix_d2eig_isppol, (2,3,ddb%mpert,3,ddb%mpert,mband,ddb%nkpt))

  ! Loop over spin index
  do isppol=1,ddb%nsppol

    bandshift = (isppol - 1)  * mband

    do iband=1,mband
      jband = bandshift + iband

      do ikpt=1,ddb%nkpt
        do ipert1=1,ddb%mpert
          do idir1=1,3
            do ipert2=1,ddb%mpert
              do idir2=1,3
                do ii=1,2
                  matrix_d2eig_isppol(ii,idir2,ipert2,idir1,ipert1,iband,ikpt)=&
                         matrix_d2eig(ii,idir2,ipert2,idir1,ipert1,jband,ikpt)
                end do
              end do
            end do
          end do
        end do
      end do


    end do ! iband

    ncerr = nf90_put_var(ncid_d2eig, nctk_idname(ncid_d2eig,&
                              'matrix_values'),&
                              matrix_d2eig_isppol,&
                              start=[1,1,1,1,1,1,1,isppol,iblok])
                              !count=[2,3,ddb%mpert,3,ddb%mpert,ddb%nband,ddb%nkpt,1,1])
                              !count=[2,3,ddb%mpert,3,ddb%mpert,ddb%nkpt,ddb%nband,1,1])
    NCF_CHECK(ncerr)

  end do  ! isppol

  ncerr = nf90_put_var(ncid_d2eig, nctk_idname(ncid_d2eig,&
                            'matrix_mask'),&
                            flg_d2eig,&
                            start=[1,1,1,1,iblok])
  NCF_CHECK(ncerr)

  ABI_SFREE(matrix_d2eig_isppol)
  ABI_SFREE(matrix_d2eig)
  ABI_SFREE(flg_d2eig)

end subroutine ddb_write_d2eig_nc
!!***

!----------------------------------------------------------------------


!!****f* m_ddb/ddb_write_d2eig_txt
!! NAME
!! ddb_write_d2eig_txt
!!
!! FUNCTION
!!  Write the eig2d data as the next block in text file format.
!!
!! INPUTS
!!  ddb_hdr=ddb header object.
!!  unddb=unit of the open ddb file in text format.
!!
!!
!! SOURCE

subroutine ddb_write_d2eig_txt(ddb, unddb, iblok)
!Arguments -------------------------------
 class(ddb_type),intent(in) :: ddb
 integer,intent(in) :: unddb
 integer,intent(in) :: iblok

!Local variables -------------------------
!scalars
 integer,parameter :: iblok_eig2d=1
 integer,parameter :: choice=2

! ************************************************************************

  ! GA: This routine is redundant with outbsd.
  !     The present implementation should replace outbsd.

 call ddb%write_block_txt(iblok,choice,ddb%nband,ddb%mpert,ddb%msize,ddb%nkpt,unddb,&
                      ddb%eig2dval(:,:,:,:), ddb%kpt(:,:))

end subroutine ddb_write_d2eig_txt
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_write_nc
!! NAME
!! ddb_write_nc
!!
!! FUNCTION
!!  Write the netndf file DDB.nc in the format of version 20230219.
!!
!! INPUTS
!!  ddb_hdr=ddb header object with open file.
!!  filename=DDB filename.
!!  comm=MPI communicator.
!!  with_psps
!!      1-> include information on pseudopoentials
!!      0-> do not include information on pseudopoentials
!!
!! SOURCE

subroutine ddb_write_nc(ddb, ddb_hdr, filename, comm, with_psps)

!Arguments -------------------------------
 class(ddb_type),intent(inout) :: ddb
 type(ddb_hdr_type),intent(inout) :: ddb_hdr
 character(len=*),intent(in) :: filename
 integer,intent(in),optional :: comm
 integer,intent(in),optional :: with_psps

!Local variables -------------------------
!scalars
 integer,parameter :: master=0
 integer :: ncid, ncerr, ncid_d0E, ncid_d1E, ncid_d2E, ncid_d3E, ncid_d2eig
 integer :: ii,iblok,iblok_d0E,iblok_d1E,iblok_d2E,iblok_d3E,iblok_d2eig
!arrays
 integer,allocatable :: flg_d1E(:,:)
 integer,allocatable :: flg_d2E(:,:,:,:)
 integer,allocatable :: flg_d3E(:,:,:,:,:,:)
 real(dp) :: qpt(3), qpts(3,3), nrms(3)
 real(dp),allocatable :: matrix_d1E(:,:,:)
 real(dp),allocatable :: matrix_d2E(:,:,:,:,:)
 real(dp),allocatable :: matrix_d3E(:,:,:,:,:,:,:)

! ************************************************************************

 if (present(comm)) then
   if (xmpi_comm_rank(comm) /= master) return
 end if

#ifdef HAVE_NETCDF

 ! =====================
 ! Header and dimensions
 ! =====================
 ! Copy types and dimensions into the header
 ddb_hdr%mpert = ddb%mpert
 call ddb_hdr%set_typ(ddb%nblok, ddb%typ)

 call ddb_hdr%open_write_nc(filename, with_psps=with_psps)
 ncid = ddb_hdr%ncid

 ! Get all group id
 ncid_d0E = nctk_idgroup(ncid, 'd0E')
 ncid_d1E = nctk_idgroup(ncid, 'd1E')
 ncid_d2E = nctk_idgroup(ncid, 'd2E')
 ncid_d3E = nctk_idgroup(ncid, 'd3E')
 ncid_d2eig = nctk_idgroup(ncid, 'd2eig')

 ! =============================
 ! Loop over block to be written
 ! =============================

 iblok_d0E = 0; iblok_d1E = 0; iblok_d2E = 0; iblok_d3E = 0; iblok_d2eig = 0

 do iblok=1,ddb%nblok

   ! ------------------------
   ! Zeroth-order derivatives
   ! ------------------------
   if (is_type_d0E(ddb%typ(iblok))) then

     iblok_d0E = iblok_d0E + 1

     ncerr = nf90_put_var(ncid_d0E, nctk_idname(ncid_d0E,&
                            'matrix_values'),&
                            ddb%val(1,1,iblok),&
                            start=[iblok_d0E])
     NCF_CHECK(ncerr)

     ncerr = nf90_put_var(ncid_d0E, nctk_idname(ncid_d0E,&
                            'matrix_mask'),&
                            ddb%flg(1,iblok),&
                            start=[iblok_d0E])
     NCF_CHECK(ncerr)

   ! -----------------------
   ! First-order derivatives
   ! -----------------------
   else if (is_type_d1E(ddb%typ(iblok))) then

     iblok_d1E = iblok_d1E + 1

     call ddb%get_d1matr(iblok, matrix_d1E, flg_d1E)

     ncerr = nf90_put_var(ncid_d1E, nctk_idname(ncid_d1E,&
                            'matrix_values'),&
                            matrix_d1E,&
                            start=[1,1,1,iblok_d1E])
     NCF_CHECK(ncerr)

     ncerr = nf90_put_var(ncid_d1E, nctk_idname(ncid_d1E,&
                            'matrix_mask'),&
                            flg_d1E,&
                            start=[1,1,iblok_d1E])
     NCF_CHECK(ncerr)
     ABI_SFREE(matrix_d1E)
     ABI_SFREE(flg_d1E)

   ! ------------------------
   ! Second-order derivatives
   ! ------------------------
   else if (is_type_d2E(ddb%typ(iblok))) then

     iblok_d2E = iblok_d2E + 1

     do ii=1,3
       qpt(ii) = ddb%qpt(ii,iblok)
     end do
     ncerr = nf90_put_var(ncid_d2E, nctk_idname(ncid_d2E,&
                            'reduced_coordinates_of_qpoints'),&
                            qpt,&
                            start=[1,iblok_d2E])
     NCF_CHECK(ncerr)

     ncerr = nf90_put_var(ncid_d2E, nctk_idname(ncid_d2E,&
                            'qpoints_normalization'),&
                            (ddb%nrm(1,iblok)),&
                            start=[iblok_d2E])
     NCF_CHECK(ncerr)

     call ddb%get_d2matr(iblok, matrix_d2E, flg_d2E)

     ncerr = nf90_put_var(ncid_d2E, nctk_idname(ncid_d2E,&
                            'matrix_values'),&
                            matrix_d2E,&
                            start=[1,1,1,1,1,iblok_d2E])
     NCF_CHECK(ncerr)

     ncerr = nf90_put_var(ncid_d2E, nctk_idname(ncid_d2E,&
                            'matrix_mask'),&
                            flg_d2E,&
                            start=[1,1,1,1,iblok_d2E])
     NCF_CHECK(ncerr)
     ABI_SFREE(matrix_d2E)
     ABI_SFREE(flg_d2E)

   ! -----------------------
   ! Third-order derivatives
   ! -----------------------
   else if (is_type_d3E(ddb%typ(iblok))) then

     iblok_d3E = iblok_d3E + 1

     do ii=1,3
       nrms(ii) = ddb%nrm(ii,iblok)
       qpts(1,ii) = ddb%qpt(ii,iblok)
       qpts(2,ii) = ddb%qpt(ii+3,iblok)
       qpts(3,ii) = ddb%qpt(ii+6,iblok)
     end do

     ncerr = nf90_put_var(ncid_d3E, nctk_idname(ncid_d3E,&
                            'reduced_coordinates_of_qpoints'),&
                            qpts,&
                            start=[1,1,iblok_d3E])
     NCF_CHECK(ncerr)

     ncerr = nf90_put_var(ncid_d3E, nctk_idname(ncid_d3E,&
                            'qpoints_normalization'),&
                            nrms,&
                            start=[1,iblok_d3E])
     NCF_CHECK(ncerr)

     call ddb%get_d3matr(iblok, matrix_d3E, flg_d3E)

     ncerr = nf90_put_var(ncid_d3E, nctk_idname(ncid_d3E,&
                            'matrix_values'),&
                            matrix_d3E,&
                            start=[1,1,1,1,1,1,1,iblok_d3E])
     NCF_CHECK(ncerr)

     ncerr = nf90_put_var(ncid_d3E, nctk_idname(ncid_d3E,&
                            'matrix_mask'),&
                            flg_d3E,&
                            start=[1,1,1,1,1,1,iblok_d3E])
     NCF_CHECK(ncerr)

     ABI_SFREE(matrix_d3E)
     ABI_SFREE(flg_d3E)

   ! ---------------------------------------
   ! Second-order derivatives of eigenvalues
   ! ---------------------------------------
   else if (is_type_d2eig(ddb%typ(iblok))) then

     iblok_d2eig = iblok_d2eig + 1

     call ddb%write_d2eig_nc(ncid_d2eig, iblok_d2eig)

   end if
 end do

#else
 ABI_ERROR("NETCDF support required to write DDB.nc file.")
#endif

end subroutine ddb_write_nc
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_read_d0E_nc
!! NAME
!! ddb_read_d0E_nc
!!
!! FUNCTION
!!  Read a DDB block containing 0th order derivatives of energy.
!!  
!!
!! INPUTS
!!  ncid=netcdf identifier of a file open in reading mode.
!!  iblok=index of the block we are setting.
!!  iblok_d0E=index of the block we are reading in the d0E group.
!!
!! OUTPUT
!!
!! SOURCE

subroutine ddb_read_d0E_nc(ddb, ncid, iblok, iblok_d0E)

!Arguments -------------------------------
!scalars
 class(ddb_type),intent(inout) :: ddb
 integer,intent(in) :: ncid,iblok,iblok_d0E

!Local variables -------------------------
!scalars
 integer :: ncid_d0E
 integer :: ncerr
!arrays
 integer :: flg(1)
 real(dp) :: val(1)

! ************************************************************************

 ncid_d0E = nctk_idgroup(ncid, 'd0E')

 ! Allocate temporary arrays
 ncerr = nf90_get_var(ncid_d0E, nctk_idname(ncid_d0E, 'matrix_values'), val, start=[1,1,iblok_d0E], count=[1,1,1])
 NCF_CHECK(ncerr)
 ddb%val(1,1,iblok) = val(1)
 ddb%val(2,1,iblok) = zero
 ncerr = nf90_get_var(ncid_d0E, nctk_idname(ncid_d0E, 'matrix_mask'), flg, start=[1,iblok_d0E], count=[1,1])
 NCF_CHECK(ncerr)
 ddb%flg(1,iblok) = flg(1)

end subroutine ddb_read_d0E_nc
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_read_d1E_nc
!! NAME
!! ddb_read_d1E_nc
!!
!! FUNCTION
!!  Read a DDB block containing 1st order derivatives of energy.
!!  
!!
!! INPUTS
!!  ncid=netcdf identifier of a file open in reading mode.
!!  iblok=index of the block we are setting.
!!  iblok_d1E=index of the block we are reading in the d1E group.
!!
!! OUTPUT
!!
!! SOURCE

subroutine ddb_read_d1E_nc(ddb, ncid, iblok, iblok_d1E)

!Arguments -------------------------------
!scalars
 class(ddb_type),intent(inout) :: ddb
 integer,intent(in) :: ncid,iblok,iblok_d1E

!Local variables -------------------------
!scalars
 integer :: ncid_d1E
 integer :: ncerr
!arrays
 integer,allocatable :: flg_d1E(:,:)
 real(dp),allocatable :: matrix_d1E(:,:,:)

! ************************************************************************

 ncid_d1E = nctk_idgroup(ncid, 'd1E')

 ! Allocate temporary arrays
 ABI_MALLOC(matrix_d1E, (2,3,ddb%mpert))
 ABI_MALLOC(flg_d1E, (3,ddb%mpert))

 ncerr = nf90_get_var(ncid_d1E, nctk_idname(ncid_d1E, 'matrix_values'), matrix_d1E, start=[1,1,1,iblok_d1E])
 NCF_CHECK(ncerr)
 ncerr = nf90_get_var(ncid_d1E, nctk_idname(ncid_d1E, 'matrix_mask'), flg_d1E, start=[1,1,iblok_d1E])
 NCF_CHECK(ncerr)

 ! Reshape
 call ddb%set_d1matr(iblok, matrix_d1E, flg_d1E)

 ! Free memory
 ABI_FREE(matrix_d1E)
 ABI_FREE(flg_d1E)

end subroutine ddb_read_d1E_nc
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_read_d2E_nc
!! NAME
!! ddb_read_d2E_nc
!!
!! FUNCTION
!!  Read a DDB block containing 2nd order derivatives of energy.
!!
!! INPUTS
!!  ncid=netcdf identifier of a file open in reading mode.
!!  iblok=index of the block we are setting.
!!  iblok_d2E=index of the block we are reading in the d2E group.
!!
!! OUTPUT
!!
!! SOURCE

subroutine ddb_read_d2E_nc(ddb, ncid, iblok, iblok_d2E)

!Arguments -------------------------------
!scalars
 class(ddb_type),intent(inout) :: ddb
 integer,intent(in) :: ncid,iblok,iblok_d2E

!Local variables -------------------------
!scalars
 integer :: ncid_d2E
 integer :: ncerr
!arrays
 real(dp) :: qpt(3)
 integer,allocatable :: flg_d2E(:,:,:,:)
 real(dp),allocatable :: matrix_d2E(:,:,:,:,:)

! ************************************************************************

 ncid_d2E = nctk_idgroup(ncid, 'd2E')

 ! Allocate temporary arrays
 ABI_MALLOC(matrix_d2E, (2,3,ddb%mpert,3,ddb%mpert))
 ABI_MALLOC(flg_d2E, (3,ddb%mpert,3,ddb%mpert))

 ncerr = nf90_get_var(ncid_d2E, nctk_idname(ncid_d2E, 'reduced_coordinates_of_qpoints'), qpt, start=[1,iblok_d2E])
 NCF_CHECK(ncerr)
 ddb%qpt(1:3,iblok) = qpt(:)
 ncerr = nf90_get_var(ncid_d2E, nctk_idname(ncid_d2E, 'qpoints_normalization'), ddb%nrm(1,iblok), start=[iblok_d2E])
 NCF_CHECK(ncerr)

 ncerr = nf90_get_var(ncid_d2E, nctk_idname(ncid_d2E, 'matrix_values'), matrix_d2E, start=[1,1,1,1,1,iblok_d2E])
 NCF_CHECK(ncerr)
 ncerr = nf90_get_var(ncid_d2E, nctk_idname(ncid_d2E, 'matrix_mask'), flg_d2E, start=[1,1,1,1,iblok_d2E])
 NCF_CHECK(ncerr)

 ! Reshape
 call ddb%set_d2matr(iblok, matrix_d2E, flg_d2E)

 ! Free memory
 ABI_FREE(matrix_d2E)
 ABI_FREE(flg_d2E)

end subroutine ddb_read_d2E_nc
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_read_d3E_nc
!! NAME
!! ddb_read_d3E_nc
!!
!! FUNCTION
!!  Read a DDB block containing 3rd order derivatives of energy.
!!
!! INPUTS
!!  ncid=netcdf identifier of a file open in reading mode.
!!  iblok=index of the block we are setting.
!!  iblok_d3E=index of the block we are reading in the d3E group.
!!
!! OUTPUT
!!
!! SOURCE

subroutine ddb_read_d3E_nc(ddb, ncid, iblok, iblok_d3E)

!Arguments -------------------------------
!scalars
 class(ddb_type),intent(inout) :: ddb
 integer,intent(in) :: ncid,iblok,iblok_d3E

!Local variables -------------------------
!scalars
 integer :: blktyp
 integer :: ncid_d3E
 integer :: ncerr
!arrays
 real(dp) :: qpt(3), nrm(3)
 real(dp),allocatable :: matrix_d3E(:,:,:,:,:,:,:)
 integer,allocatable :: flg_d3E(:,:,:,:,:,:)

! ************************************************************************

 ncid_d3E = nctk_idgroup(ncid, 'd3E')

 ! Allocate temporary arrays
 ABI_MALLOC(matrix_d3E, (2,3,ddb%mpert,3,ddb%mpert,3,ddb%mpert))
 ABI_MALLOC(flg_d3E, (3,ddb%mpert,3,ddb%mpert,3,ddb%mpert))

 ncerr = nf90_get_var(ncid_d3E, nctk_idname(ncid_d3E, 'reduced_coordinates_of_qpoints'),qpt,start=[1,1,iblok_d3E],count=[1,3,1])
 NCF_CHECK(ncerr)
 ddb%qpt(1:3,iblok) = qpt(:)

 ncerr = nf90_get_var(ncid_d3E, nctk_idname(ncid_d3E, 'reduced_coordinates_of_qpoints'),qpt,start=[2,1,iblok_d3E],count=[1,3,1])
 NCF_CHECK(ncerr)
 ddb%qpt(4:6,iblok) = qpt(:)

 ncerr = nf90_get_var(ncid_d3E, nctk_idname(ncid_d3E, 'reduced_coordinates_of_qpoints'),qpt,start=[3,1,iblok_d3E],count=[1,3,1])
 NCF_CHECK(ncerr)
 ddb%qpt(7:9,iblok) = qpt(:)

 ncerr = nf90_get_var(ncid_d3E, nctk_idname(ncid_d3E, 'qpoints_normalization'), nrm, start=[1,iblok_d3E],count=[3,1])
 NCF_CHECK(ncerr)
 ddb%nrm(:,iblok) = nrm(:)

 NCF_CHECK(nf90_get_var(ncid_d3E, nctk_idname(ncid_d3E, 'matrix_values'), matrix_d3E, start=[1,1,1,1,1,1,1,iblok_d3E]))
 NCF_CHECK(nf90_get_var(ncid_d3E, nctk_idname(ncid_d3E, 'matrix_mask'), flg_d3E, start=[1,1,1,1,1,1,iblok_d3E]))

 
 blktyp = ddb%typ(iblok) ! Save block type so it doesnt get overwritten.

 call ddb%set_d3matr(iblok, matrix_d3E, flg_d3E)

 ddb%typ(iblok) = blktyp

 ! Free memory
 ABI_FREE(matrix_d3E)
 ABI_FREE(flg_d3E)

end subroutine ddb_read_d3E_nc
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_read_d2eig_nc
!! NAME
!! ddb_read_d2eig_nc
!!
!! FUNCTION
!!  Read a DDB block containing 2nd order derivatives of eigenvalues.
!!
!! INPUTS
!!  ncid=netcdf identifier of a file open in reading mode.
!!  iblok=index of the block we are setting.
!!  iblok_d2eig=index of the block we are reading in the d2eig group.
!!
!! OUTPUT
!!
!! SOURCE

subroutine ddb_read_d2eig_nc(ddb, ncid, iblok, iblok_d2eig)

!Arguments -------------------------------
!scalars
 class(ddb_type),intent(inout) :: ddb
 integer,intent(in) :: ncid,iblok
 integer,intent(in),optional :: iblok_d2eig

!Local variables -------------------------
!scalars
 integer :: ncid_d2eig,ncerr
 integer :: nkpt_file
 integer :: nblok_d2eig
 integer :: iblok_, iblok_d2eig_
 integer :: iband, jband, bandshift, isppol, mband
 integer :: ikpt,ipert1,idir1,ipert2,idir2,ii
 character(len=500) :: msg
!arrays
 integer,allocatable :: flg_d2eig(:,:,:,:)
 real(dp) :: qpt(3)
 real(dp),allocatable :: nrm(:)
 real(dp),allocatable :: matrix_d2eig(:,:,:,:,:,:,:)
 real(dp),allocatable :: matrix_d2eig_isppol(:,:,:,:,:,:,:)

! ************************************************************************

 ncid_d2eig = nctk_idgroup(ncid, 'd2eig')

 if (present(iblok_d2eig)) then
   iblok_d2eig_= iblok_d2eig
 else
   ! Recount the blok index
   iblok_d2eig_ = 0
   do iblok_=1,iblok
     if (is_type_d2eig(ddb%typ(iblok_))) then
       iblok_d2eig_ = iblok_d2eig_ + 1
     end if
   end do
 end if

 ! Sanity check on dimensions
 if (MOD(ddb%nband, ddb%nsppol)/=0) then
    write(msg,'(a,i5,a,i5)') 'ddb was allocated with nband=',ddb%nband,&
                             ' but nsppol=',ddb%nsppol
    ABI_ERROR(msg)
 end if

 ! Read kpoints
 NCF_CHECK(nctk_get_dim(ncid, "number_of_kpoints", nkpt_file))
 ncerr = nf90_get_var(ncid, nctk_idname(ncid, 'reduced_coordinates_of_kpoints'), ddb%kpt, count=[3,nkpt_file])
 NCF_CHECK(ncerr)

 mband = ddb%nband / ddb%nsppol

 ncerr = nf90_get_var(ncid_d2eig, nctk_idname(ncid_d2eig,&
                            'reduced_coordinates_of_qpoints'),&
                            qpt,&
                            start=[1,iblok_d2eig_])
 NCF_CHECK(ncerr)

 ddb%qpt(:,iblok) = zero
 ddb%qpt(1:3,iblok) = qpt

 !ncerr = nf90_get_var(ncid_d2eig, nctk_idname(ncid_d2eig,&
 !                         'qpoints_normalization'),&
 !                         nrm(1,iblok),&
 !                         start=[iblok_d2eig_])
 NCF_CHECK(nctk_get_dim(ncid_d2eig, "number_of_d2eig_blocks", nblok_d2eig))
 ABI_MALLOC(nrm, (nblok_d2eig))
 ncerr = nf90_get_var(ncid_d2eig, nctk_idname(ncid_d2eig,&
                          'qpoints_normalization'),&
                          nrm)
 NCF_CHECK(ncerr)
 ddb%nrm(:,iblok) = zero
 ddb%nrm(1,iblok) = nrm(iblok_d2eig_)
 ABI_FREE(nrm)


 ABI_MALLOC(matrix_d2eig, (2,3,ddb%mpert,3,ddb%mpert,ddb%nband,ddb%nkpt))
 ABI_MALLOC(matrix_d2eig_isppol, (2,3,ddb%mpert,3,ddb%mpert,mband,ddb%nkpt))
 ABI_MALLOC(flg_d2eig, (3,ddb%mpert,3,ddb%mpert))

 do isppol=1,ddb%nsppol

   ncerr = nf90_get_var(ncid_d2eig, nctk_idname(ncid_d2eig,&
                            'matrix_values'),&
                            matrix_d2eig_isppol,&
                            start=[1,1,1,1,1,1,1,isppol,iblok_d2eig_])
   NCF_CHECK(ncerr)

   if (ddb%nsppol==1) then
     matrix_d2eig(:,:,:,:,:,:,:) = matrix_d2eig_isppol(:,:,:,:,:,:,:)
   else
     bandshift = (isppol - 1) * mband
     do ikpt=1,ddb%nkpt
       do iband=1,ddb%nband
         jband = bandshift + iband
         do ipert1=1,ddb%mpert
           do idir1=1,3
             do ipert2=1,ddb%mpert
               do idir2=1,3
                 do ii=1,2
                         matrix_d2eig(ii,idir2,ipert2,idir1,ipert1,jband,ikpt)=&
                  matrix_d2eig_isppol(ii,idir2,ipert2,idir1,ipert1,iband,ikpt)
                 end do
               end do
             end do
           end do
         end do
       end do
     end do
   end if
 end do

 ncerr = nf90_get_var(ncid_d2eig, nctk_idname(ncid_d2eig,&
                            'matrix_mask'),&
                            flg_d2eig,&
                            start=[1,1,1,1,iblok_d2eig_])
 NCF_CHECK(ncerr)

 ! Store values
 call ddb%set_d2eig(iblok, matrix_d2eig, flg_d2eig)

 ! Free memory
 ABI_FREE(matrix_d2eig)
 ABI_FREE(matrix_d2eig_isppol)
 ABI_FREE(flg_d2eig)

end subroutine ddb_read_d2eig_nc
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_set_d3matr
!! NAME
!! ddb_set_d3matr
!!
!! FUNCTION
!!  Set values for the third-order derivative matrix.
!!
!! INPUTS
!!  iblok=index of the block we are setting.
!!  d2matr=the third-order derivative matrix.
!!  flg=flag to indicate presence of a given element.
!!  lw=whether this type of perturbation correspond to longwave derivatives
!!
!! OUTPUT
!!
!! SOURCE

subroutine ddb_set_d3matr(ddb, iblok, d3matr, flg, lw)

!Arguments -------------------------------
!array
 class(ddb_type),intent(inout) :: ddb
 real(dp),intent(in) :: d3matr(2,3,ddb%mpert,3,ddb%mpert,3,ddb%mpert)
 integer,intent(in) :: flg(3,ddb%mpert,3,ddb%mpert,3,ddb%mpert)
!scalars
 integer,intent(in) :: iblok
 logical,intent(in),optional :: lw

!Local variables -------------------------
!scalars
 integer :: idir1,idir2,idir3,ipert1,ipert2,ipert3,index,mpert

! ************************************************************************

 mpert = ddb%mpert

 ! GA: Should be consistent among all ddb_set_dxmatr routines
 !     and alway have options to specify block type...
 ddb%typ(iblok) = BLKTYP_d3E_xx
 if (present(lw)) then
   if (lw) then
     ddb%typ(iblok) = BLKTYP_d3E_lw
   end if
 end if

 do ipert1=1,mpert
   do idir1=1,3
     do ipert2=1,mpert
       do idir2=1,3
         do ipert3=1,mpert
           do idir3=1,3

             index = idir1 + 3*((ipert1-1)+mpert*((idir2-1) + &
             & 3*((ipert2-1)+mpert*((idir3-1) + 3*(ipert3-1)))))

             ddb%flg(index,iblok) = flg(idir1,ipert1,idir2,ipert2,idir3,ipert3)
             ddb%val(:,index,iblok)= d3matr(:,idir1,ipert1,idir2,ipert2,idir3,ipert3)

           end do
         end do
       end do
     end do
   end do
 end do


end subroutine ddb_set_d3matr
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_get_d3matr
!! NAME
!! ddb_get_d3matr
!!
!! FUNCTION
!!  Transform the third-order derivative matrix
!!  from flat indices to real tensor d3matr(cplex,ncart,natom,ncart,natom,ncart,natom)
!!
!! INPUTS
!!  iblok=index of the block to get.
!!
!! OUTPUT
!!  d3matr=the third-order derivative matrix.
!!  flg=flag to indicate presence of a given element.
!!
!! SOURCE

subroutine ddb_get_d3matr(ddb, iblok, d3matr, flg)

!Arguments -------------------------------
!array
 class(ddb_type),intent(inout) :: ddb
 integer,intent(in) :: iblok
 real(dp), allocatable, intent(out) :: d3matr(:,:,:,:,:,:,:)
 integer, allocatable, intent(out) :: flg(:,:,:,:,:,:)
!scalars

!Local variables -------------------------
!scalars
 integer :: ii,idir1,idir2,idir3,ipert1,ipert2,ipert3

! ************************************************************************

 ABI_MALLOC(d3matr, (2,3,ddb%mpert,3,ddb%mpert,3,ddb%mpert))
 ABI_MALLOC(flg, (3,ddb%mpert,3,ddb%mpert,3,ddb%mpert))

 d3matr = zero

 ii=0
 do ipert3=1,ddb%mpert
   do idir3=1,3
     do ipert2=1,ddb%mpert
       do idir2=1,3
         do ipert1=1,ddb%mpert
           do idir1=1,3
             ii=ii+1
             flg(idir1,ipert1,idir2,ipert2,idir3,ipert3) = ddb%flg(ii,iblok)
             if (ddb%flg(ii,iblok) > 0) then
               d3matr(1,idir1,ipert1,idir2,ipert2,idir3,ipert3) = ddb%val(1,ii,iblok)
               d3matr(2,idir1,ipert1,idir2,ipert2,idir3,ipert3) = ddb%val(2,ii,iblok)
             end if
           end do
         end do
       end do
     end do
   end do
 end do

end subroutine ddb_get_d3matr
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_get_d2eig
!! NAME
!! ddb_get_d2eig
!!
!! FUNCTION
!!  Transform the second-order derivative matrix of eigenvalues
!!  from flat indices to real tensor d2eig(cplex,ncart,natom,ncart,natom,nband,nkpt)
!!
!! INPUTS
!!  iblok=index of the block to get.
!!
!! OUTPUTS
!!  d2eig(cplex,ncart,natom,ncart,natom,nband,nkpt) =
!!      the second-order derivative matrix of eigenvalues
!!      with the isppol index wrapped in nband.
!!  flg(ncart,natom,ncart,natom)=flag to indicate presence of a given element.
!!
!! SOURCE

subroutine ddb_get_d2eig(ddb, d2eig, flg, iblok)

!Arguments -------------------------------
!array
 class(ddb_type),intent(inout) :: ddb
 real(dp), allocatable, intent(out) :: d2eig(:,:,:,:,:,:,:)
 integer,intent(in) :: iblok
 integer, allocatable, intent(out) :: flg(:,:,:,:)
!scalars

!Local variables -------------------------
!scalars
 integer :: ii,idir1,idir2,ipert1,ipert2,iband,ikpt

! ************************************************************************

 ABI_MALLOC(d2eig, (2,3,ddb%mpert,3,ddb%mpert,ddb%nband,ddb%nkpt))
 ABI_MALLOC(flg, (3,ddb%mpert,3,ddb%mpert))

 d2eig = zero

 do ikpt=1,ddb%nkpt
   do iband=1,ddb%nband
     ii=0
     do ipert2=1,ddb%mpert
       do idir2=1,3
         do ipert1=1,ddb%mpert
           do idir1=1,3
             ii=ii+1
             flg(idir1,ipert1,idir2,ipert2) = ddb%flg(ii,iblok)
             if (ddb%flg(ii,iblok) > 0) then
               d2eig(1,idir1,ipert1,idir2,ipert2,iband,ikpt) = ddb%eig2dval(1,ii,iband,ikpt)
               d2eig(2,idir1,ipert1,idir2,ipert2,iband,ikpt) = ddb%eig2dval(2,ii,iband,ikpt)
             end if
           end do
         end do
       end do
     end do
   end do
 end do

end subroutine ddb_get_d2eig
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_set_d2eig
!! NAME
!! ddb_set_d2eig
!!
!! FUNCTION
!!  Set values for the second-order derivatives of eigenvalues.
!!
!! INPUTS
!!  iblok=index of the block we are setting.
!!  d2eig=the second-order derivative of eigenvalues.
!!  flg=flag to indicate presence of a given element.
!!
!! OUTPUT
!!
!! NOTE   
!!  Does not handle spin index. Also, sometimes, d2eig is available with flat index
!! SOURCE

subroutine ddb_set_d2eig(ddb, iblok, d2eig, flg)

!Arguments -------------------------------
!array
 class(ddb_type),intent(inout) :: ddb
 real(dp),intent(in) :: d2eig(2,3,ddb%mpert,3,ddb%mpert,ddb%nband,ddb%nkpt)
 integer,intent(in) :: flg(3,ddb%mpert,3,ddb%mpert)
!scalars
 integer,intent(in) :: iblok

!Local variables -------------------------
!scalars
 integer :: idir1,idir2,ipert1,ipert2,index,iband,ikpt,ii

! ************************************************************************

 do ikpt=1,ddb%nkpt
   do iband=1,ddb%nband
     ii = 0
     do ipert2=1,ddb%mpert
       do idir2=1,3
         do ipert1=1,ddb%mpert
           do idir1=1,3
             index=idir1+3*((ipert1-1)+ddb%mpert*((idir2-1)+3*(ipert2-1)))
             ddb%eig2dval(1,index,iband,ikpt)=d2eig(1,idir1,ipert1,idir2,ipert2,iband,ikpt)
             ddb%eig2dval(2,index,iband,ikpt)=d2eig(2,idir1,ipert1,idir2,ipert2,iband,ikpt)
             ddb%flg(index,iblok)=flg(idir1,ipert1,idir2,ipert2)
           end do !idir1
         end do !pert1
       end do !idir2
     end do !pert2
   end do  !band
 end do !kpt

end subroutine ddb_set_d2eig
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_set_d2eig_reshape
!! NAME
!! ddb_set_d2eig_reshape
!!
!! FUNCTION
!!  Set values for the second-order derivatives of eigenvalues
!!  and reshape the array received.
!!
!! INPUTS
!!  iblok=index of the block we are setting.
!!  d2eig=the second-order derivative of eigenvalues.
!!  flg=flag to indicate presence of a given element.
!!  blktyp=block type 
!!   5->real part 
!!   6->imaginary part (broadening)
!!
!! OUTPUT
!!
!! SOURCE

subroutine ddb_set_d2eig_reshape(ddb, iblok, d2eig, flg, blktyp)

!Arguments -------------------------------
!array
 class(ddb_type),intent(inout) :: ddb
 real(dp),intent(in) :: d2eig(2,ddb%nband*ddb%nsppol,ddb%nkpt,3,ddb%mpert,3,ddb%mpert)
 integer,intent(in) :: flg(3,ddb%mpert,3,ddb%mpert)
!scalars
 integer,intent(in) :: iblok
 integer,intent(in),optional :: blktyp

!Local variables -------------------------
!scalars
 integer :: idir1,idir2,ipert1,ipert2,index,iband,ikpt,ii,mband

! ************************************************************************

 ddb%typ(iblok) = BLKTYP_d2eig_re
 if (present(blktyp)) ddb%typ(iblok) = blktyp

 ! Spin polarization is wrapped in band number
 mband = ddb%nband * ddb%nsppol

 do ikpt=1,ddb%nkpt
   do iband=1,mband
     ii = 0
     do ipert2=1,ddb%mpert
       do idir2=1,3
         do ipert1=1,ddb%mpert
           do idir1=1,3
             index=idir1+3*((ipert1-1)+ddb%mpert*((idir2-1)+3*(ipert2-1)))
             ddb%flg(index,iblok)=flg(idir1,ipert1,idir2,ipert2)
             if (ddb%flg(index,iblok) > 0) then
               ddb%eig2dval(1,index,iband,ikpt)=d2eig(1,iband,ikpt,idir1,ipert1,idir2,ipert2)
               ddb%eig2dval(2,index,iband,ikpt)=d2eig(2,iband,ikpt,idir1,ipert1,idir2,ipert2)
             end if
           end do !idir1
         end do !pert1
       end do !idir2
     end do !pert2
   end do  !band
 end do !kpt

end subroutine ddb_set_d2eig_reshape
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_to_dtset
!! NAME
!! ddb_to_dtset
!!
!! FUNCTION
!!   Initialize a dataset object from ddb.
!!
!! FIXME: I don't understand the goal of this routine.
!! The dtset constructed from the DDB won't be equal to the one used to generate the DDB
!!  There's only one safe way to init dtset i.e. from file by calling the parser
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE


subroutine ddb_to_dtset(comm, dtset, filename, psps)

!Arguments -------------------------------
 integer,intent(in) :: comm
 type(dataset_type),intent(inout) :: dtset
 type(pseudopotential_type),intent(inout) :: psps
 ! type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)
 character(len=*),intent(in) :: filename
 !Local variables -------------------------
 integer :: mxnimage,unddb
!integer :: ii, nn
 type(ddb_hdr_type) :: ddb_hdr

! ************************************************************************

 ABI_UNUSED(psps%usepaw)

!Set variables
 mxnimage = 1 ! Only 1 image in the DDB

! Must read natom from the DDB before being able to allocate some arrays needed for invars9
 unddb = get_unit()
 call ddb_hdr%open_read(filename,comm=comm)
 call ddb_hdr%close()
!close ddb file, just want to read the headers
 dtset%ngfft = ddb_hdr%ngfft

! Copy scalars from ddb
 dtset%natom = ddb_hdr%natom
 dtset%mband = ddb_hdr%mband
 dtset%nkpt = ddb_hdr%nkpt
 dtset%nsym = ddb_hdr%msym
 dtset%ntypat = ddb_hdr%ntypat
 dtset%nspden = ddb_hdr%nspden
 dtset%nspinor = ddb_hdr%nspinor
 dtset%nsppol = ddb_hdr%nsppol
 dtset%occopt = ddb_hdr%occopt
 dtset%usepaw = ddb_hdr%usepaw
 dtset%intxc = ddb_hdr%intxc
 dtset%ixc = ddb_hdr%ixc
 dtset%iscf = ddb_hdr%iscf
 dtset%dilatmx = ddb_hdr%dilatmx
 dtset%ecut = ddb_hdr%ecut
 dtset%ecutsm = ddb_hdr%ecutsm
 dtset%pawecutdg = ddb_hdr%pawecutdg
 dtset%kptnrm = ddb_hdr%kptnrm
 dtset%dfpt_sciss = ddb_hdr%dfpt_sciss
 dtset%tolwfr = 1.0_dp  ! dummy
 dtset%tphysel = ddb_hdr%tphysel
 dtset%tsmear = ddb_hdr%tsmear

 ! Copy arrays from ddb
 ABI_REMALLOC(dtset%acell_orig, (3,mxnimage))
 dtset%acell_orig(1:3,1) = ddb_hdr%acell(:)

 ABI_REMALLOC(dtset%rprim_orig, (3,3,mxnimage))
 dtset%rprim_orig(1:3,1:3,1) = ddb_hdr%rprim(:,:)

 ABI_REMALLOC(dtset%rprimd_orig, (3,3,mxnimage))
 dtset%rprimd_orig(:,1,1) = ddb_hdr%rprim(:,1) * dtset%acell_orig(1,1)
 dtset%rprimd_orig(:,2,1) = ddb_hdr%rprim(:,2) * dtset%acell_orig(2,1)
 dtset%rprimd_orig(:,3,1) = ddb_hdr%rprim(:,3) * dtset%acell_orig(3,1)

 ABI_REMALLOC(dtset%amu_orig,(dtset%ntypat,mxnimage))
 dtset%amu_orig(:,1) = ddb_hdr%amu(:)

 ABI_REMALLOC(dtset%typat, (dtset%natom))
 dtset%typat(:) = ddb_hdr%typat(1:ddb_hdr%matom)

 ABI_REMALLOC(dtset%spinat, (3,dtset%natom))
 dtset%spinat(:,:) = ddb_hdr%spinat(1:3,1:ddb_hdr%matom)

 ABI_REMALLOC(dtset%xred_orig, (3,dtset%natom,mxnimage))
 dtset%xred_orig(:,:,1) = ddb_hdr%xred(1:3,1:ddb_hdr%matom)

 ABI_REMALLOC(dtset%ziontypat, (dtset%ntypat))
 dtset%ziontypat(1:ddb_hdr%mtypat) = ddb_hdr%zion(1:ddb_hdr%mtypat)

 ABI_REMALLOC(dtset%znucl,(dtset%ntypat))
 dtset%znucl(:) = ddb_hdr%znucl(1:ddb_hdr%mtypat)

 ABI_REMALLOC(dtset%nband,(dtset%nkpt))
 dtset%nband(:) = ddb_hdr%nband(1:ddb_hdr%mkpt*ddb_hdr%nsppol)

 ABI_REMALLOC(dtset%symafm,(dtset%nsym))
 dtset%symafm(:) = ddb_hdr%symafm(1:ddb_hdr%msym)

 ABI_REMALLOC(dtset%symrel, (3,3,dtset%nsym))
 dtset%symrel(:,:,:) = ddb_hdr%symrel(1:3,1:3,1:ddb_hdr%msym)

 ABI_REMALLOC(dtset%tnons,(3,dtset%nsym))
 dtset%tnons(:,:) = ddb_hdr%tnons(1:3,1:ddb_hdr%msym)

 ABI_REMALLOC(dtset%kpt,(3,dtset%nkpt))
 dtset%kpt(:,:) = ddb_hdr%kpt(1:3,1:ddb_hdr%mkpt)

 ABI_REMALLOC(dtset%wtk,(dtset%nkpt))
 dtset%wtk(:) = ddb_hdr%wtk(1:ddb_hdr%mkpt)

 ! GA: I had way too much problems implementing pawtab_copy.
 !     The script check-libpaw would report all sorts of errors.
 !     Therefore, I do a cheap copy here, copying only the relevant info.
 !call pawtab_copy(pawtab, ddb_hdr%pawtab)
 ! nn=size(pawtab)
 ! if (nn.gt.0) then
 !   do ii=1,nn
 !     pawtab(ii)%basis_size =ddb_hdr%pawtab(ii)%basis_size
 !     pawtab(ii)%lmn_size =ddb_hdr%pawtab(ii)%lmn_size
 !     pawtab(ii)%lmn2_size =ddb_hdr%pawtab(ii)%lmn2_size
 !     pawtab(ii)%rpaw =ddb_hdr%pawtab(ii)%rpaw
 !     pawtab(ii)%rshp =ddb_hdr%pawtab(ii)%rshp
 !     pawtab(ii)%shape_type =ddb_hdr%pawtab(ii)%shape_type
 !    if (allocated(pawtab(ii)%dij0)) then
 !      call alloc_copy(ddb_hdr%pawtab(ii)%dij0,  pawtab(ii)%dij0)
 !    end if
 !   end do
 ! end if

 call ddb_hdr%free()

end subroutine ddb_to_dtset
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/merge_ddb
!! NAME
!! merge_ddb
!!
!! FUNCTION
!!  Read a list of ddb files and merge them into a single ddb object.
!!
!! INPUTS
!!     nddb=number of DDBs to merge
!!     filenames=names of input DDB files
!!     outfile=name of the merged DDB file to be written
!!     dscrpt=string description of the final ddb.
!!     chkopt=option for consistency checks between DDB files
!!         (0 --> do not check header consistency between files)
!!         (1 --> check header consistency between files)
!!
!! OUTPUT
!!
!! SOURCE

subroutine merge_ddb(nddb, filenames, outfile, dscrpt, chkopt)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: nddb
 integer,intent(in) :: chkopt
!arrays
 character(len=fnlen),intent(in) :: filenames(nddb)
 character(len=fnlen),intent(in) :: outfile, dscrpt

!Local variables -------------------------
!scalars
 integer,parameter :: master=0
 integer :: iddb, iddb_mkpt, iddb_psps
 integer :: dimekb, matom, mband, mblok, mkpt, nsppol
 integer :: mtypat, lmnmax, usepaw, mblktyp, msym
 integer :: msize, msize_, mpert
 integer :: nblok, iblok, iblok1, iblok2
 integer :: comm
 logical :: eig2d, can_merge
 integer,parameter :: prtvol=0
 character(len=500) :: msg
 type(ddb_type) :: ddb, ddb2
 type(ddb_hdr_type) :: ddb_hdr, ddb_hdr2
 type(crystal_t) :: crystal

! ************************************************************************

 comm = xmpi_world

! -----------------------------------------------
! Read all headers and evaluate arrays dimensions
! -----------------------------------------------
 if (xmpi_comm_rank(comm) == master) then
   call wrtout(std_out, sjoin(ch10, " merge_ddb: Reading all headers."))
 end if

 dimekb=0 ; matom=0 ; mband=0  ; mblok=0 ; mkpt=0 ; mpert=0
 msize=0  ; mtypat=0 ; lmnmax=0 ; usepaw=0 ; mblktyp=1
 iddb_mkpt = 1 ; iddb_psps = nddb
 msym=192

 eig2d = .False.
 do iddb=1,nddb

   call ddb_hdr%open_read(filenames(iddb), comm, dimonly=1)

   matom=max(matom,ddb_hdr%matom)

   ! GA: Should get mkpt from the ddb containing d2eig, if any.
   ! In facts, since I removed comparison on the k-points
   ! in ddb_hdr_compare, I should add a check to make sure
   ! k-points are consistent when merging d2eig data.
   if (ddb_hdr%mkpt > mkpt) then
     mkpt = ddb_hdr%mkpt
     iddb_mkpt = iddb
   end if
   !mkpt=max(mkpt,ddb_hdr%mkpt)
   mtypat=max(mtypat,ddb_hdr%mtypat)
   msym=max(msym,ddb_hdr%msym)
   mband=max(mband,ddb_hdr%mband)
   dimekb=max(dimekb,ddb_hdr%psps%dimekb)
   lmnmax=max(lmnmax,ddb_hdr%psps%lmnmax)
   usepaw=max(usepaw,ddb_hdr%usepaw)
   nsppol = ddb_hdr%nsppol

   ! Count the blocks
   mblok=mblok+ddb_hdr%nblok

   ! Figure out if we are merging eig2d files
   if (is_type_d2eig(ddb_hdr%mblktyp)) then
     eig2d = .True.
   end if

   ! Figure out if we are merging d3E blocks and compute msize accordingly
   mpert = max(mpert,ddb_hdr%mpert)
   msize_ = 3 * mpert * 3 * mpert
   if (is_type_d3E(ddb_hdr%mblktyp)) msize_ = msize_ * 3 * mpert
   msize = max(msize, msize_)

   if (ddb_hdr%with_psps>0 .or. ddb_hdr%psps%usepaw > 0) then
     iddb_psps = iddb
   end if

 end do

 ddb%nsppol = nsppol

 ! ---------------
 ! Allocate arrays
 ! ---------------
 if (eig2d) then
   ! GA: We need to multiply mband by nsppol (to keep array rank below 8)
   call ddb%malloc(msize, mblok, matom, mtypat, mpert, mkpt, mband * nsppol)
 else
   call ddb%malloc(msize, mblok, matom, mtypat, mpert)
 end if

 ! -------------------------------------------------------
 ! Initialize the output ddb_hdr using the first input ddb
 ! -------------------------------------------------------

 ! GA: The last ddb is usually the one that contains the most info on pseudos
 !     however, we should check them all and figure out which one has
 !     the most info.

 call ddb_hdr%free()  ! GA: why do I need this? Try to remove
 call ddb_hdr%open_read(filenames(1), comm, &
                          matom=matom,mtypat=mtypat,mband=mband,mkpt=mkpt,&
                          msym=msym,dimekb=dimekb,lmnmax=lmnmax,usepaw=usepaw)
 call ddb_hdr%close()
 ddb_hdr%mpert = mpert
 ddb_hdr%msize = msize

 ! GA: We are setting mkpt at initialization,
 !     but netcdf file declares dimension with nkpt.
 ! TODO: Should check consistency of nkpt
 !       among of all blocks containing eig2d data.

 ! ==================
 ! Read all databases
 ! ==================

 nblok = 0

 do iddb=1,nddb

   ! Open the corresponding input DDB, and read the database file information
   write(msg, '(a,a,i6)' )ch10,' read the input derivative database number',iddb
   call wrtout(std_out,msg)

   ! Note: it is necessary to specify mkpt, otherwise the comparison will crash
   call ddb_hdr2%open_read(filenames(iddb), comm, &
                          matom=matom,mtypat=mtypat,mband=mband,mkpt=mkpt,&
                          msym=msym,dimekb=dimekb,lmnmax=lmnmax,usepaw=usepaw)
   call ddb_hdr2%close()

   if (chkopt==1)then

     ! Compare the current DDB and input DDB information.
     ! In case of an inconsistency, halt the execution.
     call wrtout(std_out, ' compare the current and input DDB information')


     ! GA: Maybe the problem is that we are comparing uninitialized pawtab
     call ddb_hdr%compare(ddb_hdr2)

   else
     ! No comparison between the current DDB and input DDB information.
     call wrtout(std_out,msg)
     write(msg, '(3a)' )&
       'No comparison/check is performed for the current and input DDB information ',ch10,&
       'because argument --nostrict was passed to the command line. '
     ABI_COMMENT(msg)
   end if

   if (chkopt==1 .or. usepaw==1) then
     call ddb_hdr%copy_missing_variables(ddb_hdr2)
   end if

   ! GA: In principle, this could be done only once,
   ! but I could not managed to do that without failing test v8[07].
   if (iddb == iddb_psps) then
     call ddb_hdr%copy_psps_from(ddb_hdr2)
   end if

   call ddb_hdr2%free()

   ! Now read the whole DDB
   call ddb2%from_file(filenames(iddb), ddb_hdr2, crystal, comm, prtvol, raw=1)
   call crystal%free()
   call ddb_hdr2%free()

   ! --------------------------------------------------------------
   ! Double loop over the blocks of the last ddb and the output ddb
   ! --------------------------------------------------------------

   do iblok2=1,ddb2%nblok

     can_merge = .false.
     do iblok1=1, nblok

       can_merge = ddb%can_merge_blocks(ddb2, iblok1, iblok2)

       if (can_merge) then
         write(msg, '(a,i5,a,a)' )' merge block #',iblok2,' from file ', filenames(iddb)
         call wrtout(std_out,msg)
         iblok = iblok1  ! Merge with previous block
         exit
       end if
     end do

     if (.not. can_merge) then
       write(msg, '(a,i5,a,a)' )' add block #',iblok2,' from file ', filenames(iddb)
       call wrtout(std_out,msg)
       nblok = nblok + 1
       iblok = nblok
     end if

     call ddb%merge_blocks(ddb2, iblok, iblok2)

   end do  ! iblok2

   ! Free memory
   call ddb2%free()

 end do  ! iddb



 ddb_hdr%nblok = nblok
 ddb%nblok = nblok
 ddb_hdr%dscrpt = dscrpt

 call ddb_hdr%set_typ(ddb%nblok, ddb%typ)

 ddb_hdr%mpert = mpert  ! This is done anyway at writing

 ! Summarize the merging phase
 write(msg, '(a,i6,a)' )' Final DDB has ',nblok,' blocks.'
 call wrtout(std_out,msg)

 ! Always use format specified with output filename.
 ! GA: Might need an extra variable to enforce a different iomode
 ddb_hdr%iomode = iomode_from_fname(outfile)

 ! GA: This is because netcdf format has more info than txt
 !     and psps might be initialized even if with_psps==0
 ! Very weird that I have to do this.
 ! TODO Do not enforce with_psps=1. Change the test reference instead.
 if (ddb_hdr%iomode/=IO_MODE_ETSF) then
   if (ddb_hdr%with_psps==0) then
     ddb_hdr%psps%dimekb = 0
     ddb_hdr%psps%lmnmax = 0
     ABI_SFREE(ddb_hdr%psps%ekb)
     ABI_SFREE(ddb_hdr%psps%indlmn)
     ABI_MALLOC(ddb_hdr%psps%ekb,(ddb_hdr%psps%dimekb,ddb_hdr%mtypat))
     ABI_MALLOC(ddb_hdr%psps%indlmn,(6,ddb_hdr%psps%lmnmax,ddb_hdr%mtypat))
     ddb_hdr%psps%ekb = zero
     ddb_hdr%psps%indlmn = zero
   end if
 end if

 ! Enforce full initialization, regardless of DDB content
 ddb_hdr%with_psps=1
 ddb_hdr%with_dfpt_vars=1

 ! Write the final ddb to file.
 if (.not. eig2d) then
   call ddb%write(ddb_hdr, outfile)
 end if

 ! =================================
 ! Second derivatives of eigenvalues
 ! =================================
 if (eig2d) then

   ! GA: Here we assume that the blocks are complete wrt perturbations.
   !     No merging of blocks occurs.
   ! TODO: Implement merging of partial d2eig blocks

   ddb%kpt(:,:) = ddb_hdr%kpt(:,:)

   ! Open the output DDB and write the header
   call ddb_hdr%open_write(outfile, with_psps=1, comm=comm)

   iblok = 0
   do iddb=1,nddb

     call ddb_hdr2%open_read(filenames(iddb), comm, &
                            matom=matom,mtypat=mtypat,mband=mband,mkpt=mkpt,&
                            msym=msym,dimekb=dimekb,lmnmax=lmnmax,usepaw=usepaw)

     do iblok2=1,ddb_hdr2%nblok

       ! Handle one block at a time
       iblok = iblok + 1
       call ddb%read_d2eig(ddb_hdr2, iblok, iblok2)
       call ddb%write_d2eig(ddb_hdr, iblok)

     end do

     call ddb_hdr2%close()  ! Close the file
     call ddb_hdr2%free()   ! Free memory

   end do

   call ddb_hdr%close()

 end if

 ! -----------
 ! Free memory
 ! -----------
 call ddb_hdr%free()
 call ddb%free()

end subroutine merge_ddb
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/lwcart
!! NAME
!! lwcart
!!
!! FUNCTION
!! Transform the 3rd-order energy derivative read from the ddb file generated by a long wave
!! calculation into cartesian coordinates, and also...
!!
!! INPUTS
!!  blkflg(3,mpert,3,mpert,3,mpert)= ( 1 if the element of the 3dte
!!   has been calculated ; 0 otherwise )
!!  d3(2,3,mpert,3,mpert,3,mpert)= matrix of the 3DTE
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space(bohr^-1)
!!  mpert =maximum number of ipert
!!  natom= number of atoms
!!  rprimd(3,3)=dimensional primitive translations (bohr)
!!
!! OUTPUT
!! carflg(3,mpert,3,mpert,3,mpert)=1 if the element of d3cart has been calculated, 0 otherwise
!! d3cart(2,3,mpert,3,mpert,3,mpert)=matrix of third-order energy derivatives in cartesian coordinates
!!
!! SOURCE

subroutine lwcart(blkflg,carflg,d3,d3cart,gprimd,mpert,natom,rprimd)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: mpert,natom
!arrays
 integer,intent(in) :: blkflg(3,mpert,3,mpert,3,mpert)
 integer,intent(out) :: carflg(3,mpert,3,mpert,3,mpert)
 real(dp),intent(in) :: d3(2,3,mpert,3,mpert,3,mpert),gprimd(3,3),rprimd(3,3)
 real(dp),intent(out) :: d3cart(2,3,mpert,3,mpert,3,mpert)

!Local variables -------------------------
!scalars
 integer :: i1dir,i1pert,i2dir,i2pert,i3dir,i3pert
 integer :: ii
!arrays
 integer :: flg1(3),flg2(3)
 real(dp) :: vec1(3),vec2(3)

! *******************************************************************

!Transform to cartesian coordinates
 d3cart(:,:,:,:,:,:,:) = d3(:,:,:,:,:,:,:)
 carflg(:,:,:,:,:,:) = 0

 do i1pert = 1, mpert
   do i2pert = 1, mpert
     do i3pert = 1, mpert

       do i2dir = 1, 3
         do i3dir = 1, 3
           do ii= 1, 2
             vec1(:) = d3cart(ii,:,i1pert,i2dir,i2pert,i3dir,i3pert)
             flg1(:) = blkflg(:,i1pert,i2dir,i2pert,i3dir,i3pert)
             call cart39(flg1,flg2,gprimd,i1pert,natom,rprimd,vec1,vec2)
             d3cart(ii,:,i1pert,i2dir,i2pert,i3dir,i3pert) = vec2(:)
             carflg(:,i1pert,i2dir,i2pert,i3dir,i3pert) = flg2(:)
           end do
         end do
       end do

       do i1dir = 1, 3
         do i3dir = 1, 3
           do ii= 1, 2
             vec1(:) = d3cart(ii,i1dir,i1pert,:,i2pert,i3dir,i3pert)
             flg1(:) = blkflg(i1dir,i1pert,:,i2pert,i3dir,i3pert)
             call cart39(flg1,flg2,gprimd,i2pert,natom,rprimd,vec1,vec2)
             d3cart(ii,i1dir,i1pert,:,i2pert,i3dir,i3pert) = vec2(:)
             carflg(i1dir,i1pert,:,i2pert,i3dir,i3pert) = flg2(:)
           end do
         end do
       end do

       do i1dir = 1, 3
         do i2dir = 1, 3
           do ii= 1, 2
             vec1(:) = d3cart(ii,i1dir,i1pert,i2dir,i2pert,:,i3pert)
             flg1(:) = blkflg(i1dir,i1pert,i2dir,i2pert,:,i3pert)
             call cart39(flg1,flg2,gprimd,i3pert,natom,rprimd,vec1,vec2)
             d3cart(ii,i1dir,i1pert,i2dir,i2pert,:,i3pert) = vec2(:)
             carflg(i1dir,i1pert,i2dir,i2pert,:,i3pert) = flg2(:)
           end do
         end do
       end do

     end do
   end do
 end do

end subroutine lwcart
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/dtqdrp
!! NAME
!! dtqdrp
!!
!! FUNCTION
!! Reads the Dynamical Quadrupole or the P^(1) Tensor
!! in the Gamma Block coming from the Derivative Data Base
!! (long wave third-order derivatives).
!!
!! INPUTS
!! blkval(2,3*mpert*3*mpert*3*mpert)= matrix of third-order energies
!! ddb_version = 8 digit integer giving date. To mantain compatibility with olderDDB files.
!! lwsym  = 0 do not symmetrize the tensor wrt efield and qvec derivative
!!             |-> 1st gradient of polarization response to atomic displacement
!!        = 1 symmetrize the tensor wrt efield and qvec derivative
!!             |-> dynamic quadrupoles
!! natom= number of atoms in unit cell
!! mpert =maximum number of ipert
!!
!! OUTPUT
!! lwtens(3,3,3,natom) = Dynamical Quadrupoles or P^(1) tensor
!!
!! SOURCE

subroutine dtqdrp(blkval,ddb_version,lwsym,mpert,natom,lwtens)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: ddb_version,lwsym,mpert,natom
!arrays
 real(dp),intent(in) :: blkval(2,3*mpert*3*mpert*3*mpert)
 real(dp),intent(out) :: lwtens(3,3,3,natom)

!Local variables -------------------------
!scalars
 integer,parameter :: cvrsio8=20100401
 integer :: elfd,iatd,iatom,qvecd
 real(dp) :: fac
 logical :: iwrite
 character(len=500) :: msg
!arrays
 real(dp) :: d3cart(2,3,mpert,3,mpert,3,mpert)

! *********************************************************************

 d3cart(1,:,:,:,:,:,:) = reshape(blkval(1,:),shape = (/3,mpert,3,mpert,3,mpert/))
 d3cart(2,:,:,:,:,:,:) = reshape(blkval(2,:),shape = (/3,mpert,3,mpert,3,mpert/))

!Define a factor to apply if DDB file has been created with the old version of 
!the longwave driver.
 if (ddb_version <= cvrsio8) then
   fac=-two
 else
   fac=one
 end if

!Extraction of quadrupoles (need symmetrization wrt qvecd and elfd)
 do iatom = 1,natom
   do iatd = 1,3
     do elfd = 1,3
       do qvecd = 1,elfd-1
         if (lwsym==1) then
           lwtens(elfd,qvecd,iatd,iatom) = fac * &
         (d3cart(2,elfd,natom+2,iatd,iatom,qvecd,natom+8)+d3cart(2,qvecd,natom+2,iatd,iatom,elfd,natom+8))
           lwtens(qvecd,elfd,iatd,iatom) = lwtens(elfd,qvecd,iatd,iatom)
         else if (lwsym==0) then
           lwtens(elfd,qvecd,iatd,iatom) = fac * d3cart(2,elfd,natom+2,iatd,iatom,qvecd,natom+8)
           lwtens(qvecd,elfd,iatd,iatom) = fac * d3cart(2,qvecd,natom+2,iatd,iatom,elfd,natom+8)
         end if
       end do
       if (lwsym==1) then
         lwtens(elfd,elfd,iatd,iatom) = fac * two*d3cart(2,elfd,natom+2,iatd,iatom,elfd,natom+8)
       else if (lwsym==0) then
         lwtens(elfd,elfd,iatd,iatom) = fac * d3cart(2,elfd,natom+2,iatd,iatom,elfd,natom+8)
       end if
     end do
   end do
 end do

 iwrite = ab_out > 0

 if (iwrite) then
   if (lwsym==1) then
     write(msg,*)' atom   dir       Qxx         Qyy         Qzz         Qyz         Qxz         Qxy'
     call wrtout([ab_out,std_out],msg)
     do iatom= 1, natom
       write(msg,'(2x,i3,3x,a3,2x,6f12.6)') iatom, 'x',lwtens(1,1,1,iatom),lwtens(2,2,1,iatom),lwtens(3,3,1,iatom), &
     & lwtens(2,3,1,iatom),lwtens(1,3,1,iatom),lwtens(1,2,1,iatom)
       call wrtout([ab_out,std_out],msg)
       write(msg,'(2x,i3,3x,a3,2x,6f12.6)') iatom, 'y',lwtens(1,1,2,iatom),lwtens(2,2,2,iatom),lwtens(3,3,2,iatom), &
     & lwtens(2,3,2,iatom),lwtens(1,3,2,iatom),lwtens(1,2,2,iatom)
       call wrtout([ab_out,std_out],msg)
       write(msg,'(2x,i3,3x,a3,2x,6f12.6)') iatom, 'z',lwtens(1,1,3,iatom),lwtens(2,2,3,iatom),lwtens(3,3,3,iatom), &
     & lwtens(2,3,3,iatom),lwtens(1,3,3,iatom),lwtens(1,2,3,iatom)
       call wrtout([ab_out,std_out],msg)
     end do
   else if (lwsym==0) then
     write(msg,*) &
   & ' atom   dir       Pxx         Pyy         Pzz         Pyz         Pxz         Pxy         Pzy         Pzx         Pyx'
     call wrtout([ab_out,std_out],msg)
     do iatom= 1, natom
       write(msg,'(2x,i3,3x,a3,2x,9f12.6)') iatom, 'x',lwtens(1,1,1,iatom),lwtens(2,2,1,iatom),lwtens(3,3,1,iatom), &
     & lwtens(2,3,1,iatom),lwtens(1,3,1,iatom),lwtens(1,2,1,iatom), &
     & lwtens(3,2,1,iatom),lwtens(3,1,1,iatom),lwtens(2,1,1,iatom)
       call wrtout([ab_out,std_out],msg)
       write(msg,'(2x,i3,3x,a3,2x,9f12.6)') iatom, 'y',lwtens(1,1,2,iatom),lwtens(2,2,2,iatom),lwtens(3,3,2,iatom), &
     & lwtens(2,3,2,iatom),lwtens(1,3,2,iatom),lwtens(1,2,2,iatom), &
     & lwtens(3,2,2,iatom),lwtens(3,1,2,iatom),lwtens(2,1,2,iatom)
       call wrtout([ab_out,std_out],msg)
       write(msg,'(2x,i3,3x,a3,2x,9f12.6)') iatom, 'z',lwtens(1,1,3,iatom),lwtens(2,2,3,iatom),lwtens(3,3,3,iatom), &
     & lwtens(2,3,3,iatom),lwtens(1,3,3,iatom),lwtens(1,2,3,iatom), &
     & lwtens(3,2,3,iatom),lwtens(3,1,3,iatom),lwtens(2,1,3,iatom)
       call wrtout([ab_out,std_out],msg)
     end do
   endif
 end if

 end subroutine dtqdrp
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_lw_copy
!! NAME
!! ddb_lw_copy
!!
!! FUNCTION
!! Copy the ddb object after reading the long wave 3rd order derivatives
!! into a new ddb_lw and resizes ddb as for 2nd order derivatives
!!
!! INPUTS
!! ddb (INOUT) = ddb block datastructure
!! mpert =maximum number of ipert
!! natom= number of atoms in unit cell
!! ntypat= number of atom types
!!
!! OUTPUT
!! ddb_lw= ddb block datastructure
!!
!! SOURCE

 subroutine ddb_lw_copy(ddb,ddb_lw,mpert,natom,ntypat)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: mpert,natom,ntypat
!arrays
 type(ddb_type),intent(inout) :: ddb
 type(ddb_type),intent(out) :: ddb_lw

!Local variables -------------------------
!scalars
 integer :: ii,nblok,nsize,cnt

! *********************************************************************

 call ddb%copy(ddb_lw)
 call ddb%free()
 nsize=3*mpert*3*mpert
 nblok=ddb_lw%nblok-count(ddb_lw%typ(:)==BLKTYP_d3E_lw)
 call ddb%malloc(nsize, nblok, natom, ntypat, mpert)

 ! Copy dimensions and static variables.
 ddb%msize = nsize
 ddb%mpert = ddb_lw%mpert
 ddb%nblok = nblok
 ddb%natom = ddb_lw%natom
 ddb%ntypat = ddb_lw%ntypat
 ddb%occopt = ddb_lw%occopt
 ddb%prtvol = ddb_lw%prtvol

 ddb%rprim = ddb_lw%rprim
 ddb%gprim = ddb_lw%gprim
 ddb%acell = ddb_lw%acell

 ! Copy the allocatable arrays.
 ddb%amu(:) = ddb_lw%amu(:)
 cnt = 0
 do ii=1,ddb_lw%nblok
   if (ddb_lw%typ(ii)/=BLKTYP_d3E_lw) then
     cnt = cnt + 1
     ddb%flg(:,cnt)   = ddb_lw%flg(1:nsize,ii)
     ddb%val(:,:,cnt) = ddb_lw%val(:,1:nsize,ii)
     ddb%typ(cnt)     = ddb_lw%typ(ii)
     ddb%nrm(:,cnt)   = ddb_lw%nrm(:,ii)
     ddb%qpt(:,cnt)   = ddb_lw%qpt(:,ii)
   end if
 end do

 end subroutine ddb_lw_copy
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/symdm9
!! NAME
!! symdm9
!!
!! FUNCTION
!! Use the set of special k points calculated by the Monkhorst & Pack Technique.
!! Check if all the information for the k points are present in
!! the DDB to determine their dynamical matrices.
!! Generate the dynamical matrices of the set of k points which
!! samples homogeneously the entire Brillouin zone.
!!
!! INPUTS
!! %flg(nsize,nblok)= flag of existence for each element of the DDB
!! %nrm(1,nblok)=norm of qpt providing normalization
!! %qpt(1<ii<9,nblok)=q vector of a phonon mode (ii=1,2,3)
!! %typ(nblok)=1 or 2 depending on non-stationary or stationary block 3 for third order derivatives
!! %val(2,3*mpert*3*mpert,nblok)= all the dynamical matrices
!! gprim(3,3)=dimensionless primitive translations in reciprocal space
!! indsym = mapping of atoms under symops
!! mpert =maximum number of ipert
!! natom=number of atoms in unit cell
!! %nblok=number of blocks in the DDB
!! nqpt=number of special q points
!! nsym=number of space group symmetries
!! rfmeth =
!!   1 or -1 if non-stationary block
!!   2 or -2 if stationary block
!!   3 or -3 if third order derivatives
!!   positive if symmetries are used to set elements to zero whenever possible, negative to prevent this to happen.
!! rprim(3,3)=dimensionless primitive translations in real space
!! spqpt(3,nqpt)=set of special q points generated by the Monkhorst & Pack Method
!! symrec(3,3,nsym)=3x3 matrices of the group symmetries (reciprocal space)
!! symrel(3,3,nsym)=3x3 matrices of the group symmetries (real space)
!! comm=MPI communicator.
!!
!! OUTPUT
!! dynmat(2,3,natom,3,natom,nqpt)=dynamical matrices relative to the q points of the B.Z. sampling
!! [qmissing]=Allocatable array with the indices of the q-points in the BZ that could not be obtained
!!    by symmetry. If qmissing is present, the routine does not stop if the full BZ cannot be reconstructed.
!!    The caller is responsible for filling the missing entries.
!!
!! TODO
!!   * A full description of the inputs should be included
!!
!! NOTES
!!   Time-reversal symmetry is always assumed
!!
!! SOURCE

subroutine symdm9(ddb, dynmat, gprim, indsym, mpert, natom, nqpt, nsym, rfmeth,&
                  rprim, spqpt, symrec, symrel, comm, qmissing)

!Arguments -------------------------------
!scalars
 type(ddb_type),intent(in) :: ddb
 integer,intent(in) :: mpert,natom,nqpt,nsym,rfmeth,comm
!arrays
 !integer,intent(in) :: blkflg(3,mpert,3,mpert,nblok),blktyp(nblok)
 integer,intent(in) :: indsym(4,nsym,natom),symrec(3,3,nsym),symrel(3,3,nsym)
 integer,allocatable,optional,intent(out) :: qmissing(:)
 !real(dp),intent(in) :: blknrm(3,nblok),blkqpt(9,nblok)
 !real(dp),intent(in) :: blkval(2,3*mpert*3*mpert,nblok),gprim(3,3),rprim(3,3)
 real(dp),intent(in) :: gprim(3,3),rprim(3,3)
 real(dp),intent(in) :: spqpt(3,nqpt)
 real(dp),intent(out) :: dynmat(2,3,natom,3,natom,nqpt)

!Local variables -------------------------
!scalars
 integer :: ia,ib,iblok,idir1,idir2,ii,ipert1,ipert2,iqpt,isym,jj,kk,ll
 integer :: mu,nu,q1,q2,nqmiss,nprocs,my_rank,ierr,index
 real(dp),parameter :: tol=2.d-8
!tolerance for equality of q points between those of the DDB and those of the sampling grid
 real(dp) :: arg1,arg2,im,re,sumi,sumr
 logical :: allow_qmiss
 character(len=500) :: msg
!arrays
 integer,allocatable :: qtest(:,:)
 integer :: qmiss_(nqpt)
 real(dp) :: qq(3),qsym(6),ss(3,3)
 real(dp),allocatable :: ddd(:,:,:,:,:)

! *********************************************************************

 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 ! Initialize output (some q-points might not be reconstructed if qmissing is present)
 dynmat = zero
 allow_qmiss = (present(qmissing))

 ABI_MALLOC(ddd,(2,3,natom,3,natom))
 ! Check if the blkqpt points and their symmetrics are sufficient
 ! in the DDB to retrieve all the q points of the B.Z. sampling

 !Initialization of a test variable
 ! qtest(iqpt,1)=iblok
 ! qtest(iqpt,2)=isym
 ! qtest(iqpt,3)=time_reversal
 ABI_MALLOC(qtest,(nqpt,3))
 do iqpt=1,nqpt
   qtest(iqpt,1)=0
 end do

 !Q points coming from the DDB
 !write(std_out,*)' Nbr. of Blocks -> ',nblok
 ! TODO: This part scales badly with nblock/nqpt
 ! One could use listkk or rearrange the loop so that iqpt comes first and then MPI-parallelize.

 do iblok=1,ddb%nblok

   if (abs(ddb%typ(iblok)) == abs(rfmeth)) then
     qq(1)=ddb%qpt(1,iblok)/ddb%nrm(1,iblok)
     qq(2)=ddb%qpt(2,iblok)/ddb%nrm(1,iblok)
     qq(3)=ddb%qpt(3,iblok)/ddb%nrm(1,iblok)

     ! Calculation of the symmetric points (including Time Reversal)
     do isym=1,nsym
       qsym(1)=qq(1)*symrec(1,1,isym)+qq(2)*symrec(1,2,isym)+qq(3)*symrec(1,3,isym)
       qsym(2)=qq(1)*symrec(2,1,isym)+qq(2)*symrec(2,2,isym)+qq(3)*symrec(2,3,isym)
       qsym(3)=qq(1)*symrec(3,1,isym)+qq(2)*symrec(3,2,isym)+qq(3)*symrec(3,3,isym)

       ! Dont forget the Time Reversal symmetry
       qsym(4)=-qq(1)*symrec(1,1,isym)-qq(2)*symrec(1,2,isym)-qq(3)*symrec(1,3,isym)
       qsym(5)=-qq(1)*symrec(2,1,isym)-qq(2)*symrec(2,2,isym)-qq(3)*symrec(2,3,isym)
       qsym(6)=-qq(1)*symrec(3,1,isym)-qq(2)*symrec(3,2,isym)-qq(3)*symrec(3,3,isym)

       ! Comparison between the q points and their symmetric points
       ! and the set of q points which samples the entire Brillouin zone
       do iqpt=1,nqpt

         if (mod(abs(spqpt(1,iqpt)-qsym(1))+tol,1._dp)<2*tol)then
           if (mod(abs(spqpt(2,iqpt)-qsym(2))+tol,1._dp)<2*tol)then
             if (mod(abs(spqpt(3,iqpt)-qsym(3))+tol,1._dp)<2*tol)then

               ! write(std_out,*)' q point from the DDB ! '
               ! write(std_out,*)' block -> ',iblok
               ! write(std_out,*)' sym.  -> ',isym
               ! write(std_out,*)' No Time Reversal '
               ! write(std_out,*)'(',qsym(1),',',qsym(2),',',qsym(3),')'
               ! write(std_out,*)' '
               qtest(iqpt,1)=iblok
               qtest(iqpt,2)=isym
               qtest(iqpt,3)=0
             end if
           end if
         end if

         if (mod(abs(spqpt(1,iqpt)-qsym(4))+tol,1._dp)<2*tol)then
           if (mod(abs(spqpt(2,iqpt)-qsym(5))+tol,1._dp)<2*tol)then
             if (mod(abs(spqpt(3,iqpt)-qsym(6))+tol,1._dp)<2*tol)then

               ! write(std_out,*)' q point from the DDB ! '
               ! write(std_out,*)' block -> ',iblok
               ! write(std_out,*)' sym.  -> ',isym
               ! write(std_out,*)' Time Reversal '
               ! write(std_out,*)'(',qsym(4),',',qsym(5),',',qsym(6),')'
               ! write(std_out,*)' '

               qtest(iqpt,1)=iblok
               qtest(iqpt,2)=isym
               qtest(iqpt,3)=1
             end if
           end if
         end if

       end do ! iqpt
     end do ! isym

   end if
 end do ! iblok

! Check if all the information relatives to the q points sampling are found in the DDB;
! if not => stop message
 nqmiss = 0
 do iqpt=1,nqpt
   if (qtest(iqpt,1)==0) then
     nqmiss = nqmiss + 1
     qmiss_(nqmiss) = iqpt
     write(msg, '(3a)' )&
      ' symdm9: the bloks found in the DDB are characterized',ch10,&
      '  by the following wavevectors :'
     call wrtout(std_out,msg)
     do iblok=1,ddb%nblok
       write(msg, '(a,4d20.12)')' ',ddb%qpt(1,iblok),ddb%qpt(2,iblok),ddb%qpt(3,iblok),ddb%nrm(1,iblok)
       call wrtout(std_out,msg)
     end do
     write(msg, '(a,a,a,i0,a,a,a,3es16.6,a,a,a,a)' )&
      'Information is missing in the DDB.',ch10,&
      'The dynamical matrix number ',iqpt,' cannot be built,',ch10,&
      'since no block with qpt:',spqpt(1:3,iqpt),ch10,&
      'has been found.',ch10,&
      'Action: add the required block in the DDB, or modify the q-mesh your input file.'
     if (.not.allow_qmiss) then
       ABI_ERROR(msg)
     else
       !continue
       ABI_COMMENT(msg)
     end if
   end if
 end do

 ! Will return a list with the index of the q-points that could not be symmetrized.
 if (allow_qmiss) then
   ABI_MALLOC(qmissing, (nqmiss))
   if (nqmiss > 0) qmissing = qmiss_(1:nqmiss)
 end if

 ! Generation of the dynamical matrices relative to the q points
 ! of the set which samples the entire Brillouin zone
 do iqpt=1,nqpt
   if (mod(iqpt, nprocs) /= my_rank) cycle ! mpi-parallelism

   q1=qtest(iqpt,1)
   q2=qtest(iqpt,2)
   ! Skip this q-point if don't have enough info and allow_qmiss
   if (allow_qmiss .and. q1==0) cycle

   ! Check if the symmetry accompagnied with time reversal : q <- -q
   if (qtest(iqpt,3)==0) then
     do ii=1,3
       qq(ii)=ddb%qpt(ii,q1)/ddb%nrm(1,q1)
       do jj=1,3
         ss(ii,jj)=zero
         do kk=1,3
           do ll=1,3
             ss(ii,jj)=ss(ii,jj)+rprim(ii,kk)*gprim(jj,ll)*symrel(kk,ll,q2)
           end do
         end do
       end do
     end do
   else
     do ii=1,3
       qq(ii)=-ddb%qpt(ii,q1)/ddb%nrm(1,q1)
       do jj=1,3
         ss(ii,jj)=zero
         do kk=1,3
           do ll=1,3
             ss(ii,jj)=ss(ii,jj)+rprim(ii,kk)*gprim(jj,ll)*symrel(kk,ll,q2)
           end do
         end do
       end do
     end do
   end if

   ! Check whether all the information is contained in the DDB
   do ipert2=1,natom
     do idir2=1,3
       do ipert1=1,natom
         do idir1=1,3
           index = idir1+ 3*((ipert1-1)+ddb%mpert*((idir2-1)+3*(ipert2-1)))
           !if(ddb%flg(idir1,ipert1,idir2,ipert2,q1)/=1)then
           if(ddb%flg(index,q1)/=1)then
             write(msg, '(a,a,a,i0,a,a,a,4(i0,1x),a,a,a,a)' )&
             'Elements are missing in the DDB.',ch10,&
             'In block iq1: ',q1,' the following element is missing: ',ch10,&
             '(idir1, ipert1, idir2, ipert2): ',idir1,ipert1,idir2,ipert2,ch10,&
             'Action: add the required information in the DDB with mrgddb,',ch10,&
             'and/or check that all irreducible perturbations have been computed.'
             ABI_ERROR(msg)
           end if
         end do
       end do
     end do
   end do

   ! Read the dynamical matrices in the DDB
   do ipert2=1,natom
     do idir2=1,3
       do ipert1=1,natom
         do idir1=1,3
           ddd(:,idir1,ipert1,idir2,ipert2)=ddb%val(:,idir1+3*(ipert1-1+mpert*(idir2-1+3*(ipert2-1))),q1)
         end do
       end do
     end do
   end do

   ! Calculation of the dynamical matrix of a symmetrical q point
   do ia=1,natom
     do ib=1,natom
       ! write(std_out,*)'atom-> ',ia,indsym(4,q2,ia); write(std_out,*)'atom-> ',ib,indsym(4,q2,ib)
       arg1=two_pi*(qq(1)*indsym(1,q2,ia)+qq(2)*indsym(2,q2,ia)+qq(3)*indsym(3,q2,ia))
       arg2=two_pi*(qq(1)*indsym(1,q2,ib)+qq(2)*indsym(2,q2,ib)+qq(3)*indsym(3,q2,ib))
       re=cos(arg1)*cos(arg2)+sin(arg1)*sin(arg2)
       im=cos(arg2)*sin(arg1)-cos(arg1)*sin(arg2)
       do mu=1,3
         do nu=1,3
           sumr=zero
           sumi=zero
           do ii=1,3
             do jj=1,3
               ! If there is Time Reversal : D.M. <- Complex Conjugate D.M.
               if (qtest(iqpt,3)==0) then
                 sumr=sumr+ss(mu,ii)*ss(nu,jj)*ddd(1,ii,indsym(4,q2,ia),jj,indsym(4,q2,ib))
                 sumi=sumi+ss(mu,ii)*ss(nu,jj)*ddd(2,ii,indsym(4,q2,ia),jj,indsym(4,q2,ib))
               else
                 sumr=sumr+ss(mu,ii)*ss(nu,jj)*ddd(1,ii,indsym(4,q2,ia),jj,indsym(4,q2,ib))
                 sumi=sumi-ss(mu,ii)*ss(nu,jj)*ddd(2,ii,indsym(4,q2,ia),jj,indsym(4,q2,ib))
               end if
             end do
           end do

           ! Dynmat -> Dynamical Matrix for the q point of the sampling
           ! write(std_out,*)' Sumr -> ',mu,nu,sumr; write(std_out,*)' Sumi -> ',mu,nu,sumi
           dynmat(1,mu,ia,nu,ib,iqpt)=re*sumr-im*sumi
           dynmat(2,mu,ia,nu,ib,iqpt)=re*sumi+im*sumr
         end do ! coordinates
       end do

     end do ! ia atoms
   end do ! ib atoms
 end do ! q points of the sampling

 ABI_FREE(ddd)
 ABI_FREE(qtest)


 call xmpi_sum(dynmat, comm, ierr)

end subroutine symdm9
!!***

!----------------------------------------------------------------------

end module m_ddb
!!***

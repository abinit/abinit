!{\src2tex{textfont=tt}}
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
!! Copyright (C) 2011-2019 ABINIT group (MJV, XG, MT, MM, MVeithen, MG, PB, JCC, SP)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_ddb

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
 use m_ddb_hdr
 use m_dtset

 use m_io_tools,       only : file_exists
 use defs_datatypes,   only : pseudopotential_type
 use m_fstrings,       only : sjoin, itoa, ktoa
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

 public :: dfptnl_doutput   ! Write the matrix of third-order derivatives to the output file and the DDB

 public :: rdddb9           ! This routine reads the derivative database entirely,
 public :: nlopt            ! Output of all quantities related to third-order derivatives of the energy.
 public :: chkin9
 public :: carttransf       ! Transform a second-derivative matrix (EIG2D) from reduced
                            ! coordinates to cartesian coordinates.
#ifdef MR_DEV
 public :: dfpt_lw_doutput   ! Write the matrix of third-order derivatives to the output file and the DDB
#endif

 integer,public,parameter :: DDB_VERSION=100401
 ! DDB Version number.
 ! TODO: Remove other occurrences of this magic number.
 ! Postponed to avoid conflicts with Samuel's branch in dfpt_looppert

 real(dp),public,parameter :: DDB_QTOL=2.0d-8
 ! Tolerance for the identification of two wavevectors
!!***

!!****t* m_ddb/ddb_type
!! NAME
!! ddb_type
!!
!! FUNCTION
!!  Provides methods to extract and postoprocess the results in the derivative database (DDB)
!!
!! SOURCE

 type,public :: ddb_type

  integer :: msize
  ! Maximum size of dynamical matrices and other perturbations (ddk, dde...)

  integer :: mpert
  ! TODO: Write function that returns mpert from natom!

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

  ! These values are used to call the anaddb routines that don't use rprimd, gprimd.
  real(dp) :: rprim(3,3)
  real(dp) :: gprim(3,3)
  real(dp) :: acell(3)

  integer,allocatable :: flg(:,:)
  ! flg(msize,nblok)
  ! flag to indicate presence of a given block

  integer,allocatable :: typ(:)
  ! typ(nblok)
  ! type of each block - ddk, dde, phonon etc...

  real(dp),allocatable :: amu(:)
  ! amu(ntypat)
  ! mass of the atoms (atomic mass unit)

  real(dp),allocatable :: nrm(:,:)
  ! nrm(3,nblok)
  ! norm of the q-points for each block - can be 0 to indicate a direction of approach to gamma

  real(dp),allocatable :: qpt(:,:)
  ! qpt(9,nblok)
  ! q-point vector in reciprocal space (reduced lattice coordinates) for each block

  real(dp),allocatable :: val(:,:,:)
  ! val(2,msize,nblok)
  ! values of the second energy derivatives in each block

  contains

    procedure :: free => ddb_free
     ! Free dynamic memory.

    procedure :: malloc => ddb_malloc
     ! Allocate dynamic memory

    procedure :: bcast => ddb_bcast
     ! Broadcast the object.

    procedure :: get_etotal => ddb_get_etotal
     ! Read the GS total energy.

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

    procedure :: write_block => ddb_write_block
    !procedure :: write_blok => ddb_write_blok
     ! Writes blocks of data in the DDBs.

    !procedure :: read_blok => ddb_read_block
    procedure :: read_block => ddb_read_block
     ! This routine reads blocks of data in the DDBs.

    procedure :: get_block => ddb_get_block
    !procedure :: gtblk9 => ddb_get_block
     ! Finds the block containing the derivatives of the total energy.

 end type ddb_type

 public :: ddb_from_file            ! Construct the object from the DDB file.
 public :: ddb_copy                 ! Copy the object.
 public :: ddb_to_dtset             ! Transfer ddb_hdr to dtset datatype

 public :: mblktyp1                 ! This routine merges the derivative databases of type 0-4:

 ! TODO: This routine is deprecated and will be removed
 public :: mblktyp5                 ! This routine merges the derivative databases of type 5:

 ! TODO: Add option to change amu.
 !public :: ddb_change_amu
 !public :: ddb_print
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

!!****f* m_ddb/ddb_free
!! NAME
!! ddb_free
!!
!! FUNCTION
!!  Clean and deallocate types for the ddb_type structure
!!
!! PARENTS
!!      anaddb,ddb_interpolate,dfpt_looppert,dfptnl_doutput,eph,gstate
!!      m_effective_potential_file,m_gruneisen,mblktyp1,mblktyp5,thmeig
!!
!! CHILDREN
!!      ddb_hdr_free,ddb_hdr_open_read
!!
!! SOURCE

subroutine ddb_free(ddb)

!Arguments ------------------------------------
 class(ddb_type),intent(inout) :: ddb

! ************************************************************************

 !integer
 ABI_SFREE(ddb%flg)
 ABI_SFREE(ddb%typ)

 ! real
 ABI_SFREE(ddb%amu)
 ABI_SFREE(ddb%nrm)
 ABI_SFREE(ddb%qpt)
 ABI_SFREE(ddb%val)

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
!! PARENTS
!!
!! CHILDREN
!!      ddb_hdr_free,ddb_hdr_open_read
!!
!! SOURCE

subroutine ddb_copy(iddb, oddb)

!Arguments ------------------------------------
!array
 type(ddb_type),intent(in) :: iddb
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
!!  Allocate dynamic memory
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      ddb_interpolate,dfptnl_doutput,gstate,m_ddb,mblktyp1,mblktyp5,thmeig
!!
!! CHILDREN
!!      ddb_hdr_free,ddb_hdr_open_read
!!
!! SOURCE

subroutine ddb_malloc(ddb, msize, nblok, natom, ntypat)

!Arguments ------------------------------------
!array
 integer,intent(in) :: msize,nblok,natom,ntypat
 class(ddb_type),intent(inout) :: ddb

! ************************************************************************

 ! FIXME
 ! This is done in rdddb9 but only by the master node!
 ! Should rationalize the different steps
 ddb%msize = msize
 ddb%nblok = nblok
 ddb%natom = natom
 ddb%mpert = natom+MPERT_MAX
 ddb%ntypat = ntypat

 ! integer
 ABI_CALLOC(ddb%flg,(msize,nblok))
 ABI_CALLOC(ddb%typ,(nblok))

 ! real
 ABI_MALLOC(ddb%amu,(ntypat))
 ABI_MALLOC(ddb%nrm,(3,nblok))
 ABI_MALLOC(ddb%qpt,(9,nblok))
 ABI_MALLOC(ddb%val,(2,msize,nblok))
 ddb%val = huge(one)

end subroutine ddb_malloc
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
!!   master=Rank of Master
!!   comm=MPI communicator
!!
!! SIDE EFFECTS
!!   Ddb<type(ddb_type)>= Input if node is master, other nodes returns with a completely initialized instance.
!!
!! PARENTS
!!      m_ddb
!!
!! CHILDREN
!!      ddb_hdr_free,ddb_hdr_open_read
!!
!! SOURCE

subroutine ddb_bcast(ddb, master, comm)

!Arguments ------------------------------------
!array
 class(ddb_type),intent(inout) :: ddb
 integer, intent(in) :: master,comm

!Local variables-------------------------------
!scalars
 integer :: ierr
! *************************************************************************

 if (xmpi_comm_size(comm) == 1) return

 DBG_ENTER("COLL")

 ! Transmit dimensions and static variables.
 call xmpi_bcast(ddb%msize, master, comm, ierr)
 call xmpi_bcast(ddb%nblok, master, comm, ierr)
 call xmpi_bcast(ddb%natom, master, comm, ierr)
 call xmpi_bcast(ddb%ntypat, master, comm, ierr)

 call xmpi_bcast(ddb%occopt, master, comm, ierr)
 call xmpi_bcast(ddb%prtvol, master, comm, ierr)

 !real
 call xmpi_bcast(ddb%rprim, master, comm, ierr)
 call xmpi_bcast(ddb%gprim, master, comm, ierr)
 call xmpi_bcast(ddb%acell, master, comm, ierr)

 ! Allocate arrays on the other nodes.
 if (xmpi_comm_rank(comm) /= master) then
   call ddb%malloc(ddb%msize, ddb%nblok, ddb%natom, ddb%ntypat)
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
!! This routine (get block) finds the block that contains the
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
!! rfphon(4) = 1=> response to phonons (for the four possible derivatives.
!!             Two should be used for a second derivative of total energy)
!! rfelfd(4) = 1=> d/dk, 2=> electric field only, 3=> both (see comment on rfphon)
!! rfstrs(4) = 1=> uniaxial stresses, 2=> shear stresses, 3=> both (see comment on rfphon)
!! rftyp =
!!   0 => total energy
!!   1 => non-stationary formulation of the 2nd derivative
!!   2 => stationary formulation of the 2nd derivative
!!   3 => third derivative of total energy
!!   4 => first-order derivatives of total energy
!!
!! OUTPUT
!! iblok= number of the block that corresponds to the specifications
!!
!! PARENTS
!!      anaddb,ddb_interpolate,m_ddb,m_effective_potential_file,m_phonons
!!      thmeig
!!
!! CHILDREN
!!      ddb_hdr_free,ddb_hdr_open_read
!!
!! SOURCE


subroutine ddb_get_block(ddb,iblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,rftyp)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: rftyp
 integer,intent(out) :: iblok
 class(ddb_type),intent(in) :: ddb
!arrays
 integer,intent(in) :: rfelfd(4),rfphon(4),rfstrs(4)
 real(dp),intent(inout) :: qphnrm(3),qphon(3,3)

!Local variables -------------------------
!scalars
 integer :: blkgam,ider,idir,idir1,idir2,idir3,ii,index,ipert,ipert1,ipert2
 integer :: ipert3,nder,ok,mpert,natom
 character(len=500) :: message
!arrays
 integer :: gamma(3)
 integer,allocatable :: worki(:,:)
 real(dp) :: qpt(3)

! *********************************************************************

 mpert = ddb%mpert
 natom = ddb%natom

!Get the number of derivative
 if(rftyp==1.or.rftyp==2)then
   nder=2
 else if(rftyp==3)then
   nder=3
 else if(rftyp==0)then
   nder=0
 else if(rftyp==4)then
   nder=1
 else
   write(message, '(a,i0,a)')' rftyp is equal to ',rftyp,'. The only allowed values are 0, 1, 2, 3 or 4.'
   MSG_BUG(message)
 end if

!In case of a second-derivative, a second phonon wavevector is provided.
 if(nder==2)then
   do ii=1,3
     qphon(ii,2)=-qphon(ii,1)
   end do
   qphnrm(2)=qphnrm(1)
 end if

!In case of a third derivative, the sum of wavevectors to gamma is checked
 if (nder == 3) then
   qpt(:) = qphon(:,1)/qphnrm(1) + qphon(:,2)/qphnrm(2) + qphon(:,3)/qphnrm(3)
   call gamma9(gamma(nder),qpt,qphnrm(1),DDB_QTOL)
   if (gamma(nder) == 0) then
     write(message,'(a,a,a)')&
&     'the sum of the wavevectors of the third-order energy is ',ch10,&
&     'not equal to zero'
     MSG_ERROR(message)
   end if
 end if

!Check the validity of the requirement
 do ider=1,nder
   ! Identifies if qphon is at gamma
   call gamma9(gamma(ider),qphon(1:3,ider),qphnrm(ider),DDB_QTOL)

   if(gamma(ider)==0)then
     if(rfstrs(ider)/=0.or.rfelfd(ider)/=0)then
       write(message, '(a,a)' )&
&       'Not yet able to handle stresses or electric fields',ch10,&
&       'with non-zero wavevector.'
       MSG_BUG(message)
     end if
   end if
 end do

!Initialise the perturbation table
 ABI_MALLOC(worki,(mpert,4))
 worki(:,1:nder)=0

!Build the perturbation table
 do ider=1,nder
!  First the phonons
   if(rfphon(ider)==1)then
     do ipert=1,natom
       worki(ipert,ider)=1
     end do
   end if
!  Then the d/dk
   if(rfelfd(ider)==1.or.rfelfd(ider)==3)then
     worki(natom+1,ider)=1
   end if
!  Then the electric field
   if(rfelfd(ider)==2.or.rfelfd(ider)==3)then
     worki(natom+2,ider)=1
   end if
!  Then the uniaxial stress
   if(rfstrs(ider)==1.or.rfstrs(ider)==3)then
     worki(natom+3,ider)=1
   end if
!  At last, the shear stress
   if(rfstrs(ider)==2.or.rfstrs(ider)==3)then
     worki(natom+4,ider)=1
   end if
 end do

!Examine every blok :
 do iblok=1,ddb%nblok

!  If this variable is still 1 at the end of the examination, the blok is the good one...
   ok=1

!  Check the type
   if(rftyp/=ddb%typ(iblok)) ok=0

!  Check the wavevector
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
           if( abs( ddb%qpt(idir+3*(ider-1),iblok)/ddb%nrm(ider,iblok) - &
&           qphon(idir,ider)/qphnrm(ider) )>DDB_QTOL )then
             ok=0
           end if ! qphon
         end do ! idir
       end do ! nder
     end if  ! nder

   end if ! ok

!  Check if there is enough information in this blok
   if( ok==1 )then

     if (nder == 0) then
       if (ddb%flg(1,iblok) /= 1) then
         ok = 0
         if (ddb%prtvol > 1) then
           write(message,'(a,i0,3a)' )&
           'The block ',iblok,' does not match the requirement',ch10,&
           'because it lacks the total energy'
           MSG_COMMENT(message)
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
&                         3*((ipert1 - 1) + mpert*((idir2 - 1) + &
&                         3*((ipert2 -1 ) + mpert*((idir3 - 1) + 3*(ipert3 - 1)))))
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

!  Now that everything has been checked, eventually end the search
   if(ok==1)exit
 end do

 if(ok==0)then
   iblok=0

   if (ddb%prtvol > 1) then
     write(message, '(3a)' )&
&     ' gtblk9 : ',ch10,&
&     '  Unable to find block corresponding to the following specifications :'
     call wrtout(std_out,message,'COLL')
     write(message, '(a,i3)' )' Type (rfmeth) =',rftyp
     call wrtout(std_out,message,'COLL')
     write(message, '(a)' ) ' ider qphon(3)         qphnrm   rfphon rfelfd rfstrs'
     call wrtout(std_out,message,'COLL')
     do ider=1,nder
       write(message, '(i4,4f6.2,3i7)' )&
       ider,(qphon(ii,ider),ii=1,3),qphnrm(ider),rfphon(ider),rfelfd(ider),rfstrs(ider)
       call wrtout(std_out,message,'COLL')
     end do
   end if
 end if

 if (ok==1 .and. ddb%prtvol > 1) then
   write(message,'(a,i0,2a)')' gtblk9: found block number ',iblok,' agree with',' specifications '
   call wrtout(std_out,message,'COLL')
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
!! qphon(3)=wavevector
!! qphnrm=normalisation factor
!! qtol=tolerance
!!
!! OUTPUT
!! gamma= if 1, means that the wavevector is indeed at Gamma otherwise 0.
!!
!! PARENTS
!!      m_ddb
!!
!! CHILDREN
!!      ddb_hdr_free,ddb_hdr_open_read
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

!!****f* m_db_blk/ddb_read_block
!!
!! NAME
!! ddb_read_block
!!
!! FUNCTION
!! This routine reads blocks of data in the DDBs.
!!
!! INPUTS
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
!!  used in case of a third order derivative of total energy (three wavevector could be present)
!! ddb%nrm(3)=normalization factors for the three allowed wavevectors.
!! ddb%val(2,msize)=real(dp), complex, value of the matrix elements that are present in the data block
!! [blkval2(2,msize,mband,nkpt)]= value of the matrix elements that are present in a block of EIGR2D/EIGI2D
!!
!! NOTES
!! only executed by one processor.
!!
!! PARENTS
!!      m_ddb,mblktyp1,mblktyp5,thmeig
!!
!! CHILDREN
!!      ddb_hdr_free,ddb_hdr_open_read
!!
!! SOURCE

subroutine ddb_read_block(ddb,iblok,mband,mpert,msize,nkpt,nunit,&
&     blkval2,kpt) !optional

!Arguments -------------------------------
!scalars
 integer,intent(in) :: mband,mpert,msize,nkpt,nunit
 integer, intent(in) :: iblok
 class(ddb_type),intent(inout) :: ddb
!arrays
 real(dp),intent(out),optional :: kpt(3,nkpt)
 real(dp),intent(out),optional :: blkval2(2,msize,mband,nkpt)

!Local variables -------------------------
!scalars
 integer :: band,iband,idir1,idir2,idir3,ii,ikpt,index,ipert1,ipert2,ipert3,nelmts
 real(dp) :: ai,ar
 character(len=32) :: name
 character(len=500) :: message

! *********************************************************************

!Zero every flag
 ddb%flg(1:msize, iblok)=0
 if(present(blkval2))blkval2(:,:,:,:)=zero
 if(present(kpt))kpt(:,:)=zero

!Read the block type and number of elements
 read(nunit,*)
 read(nunit, '(a32,12x,i8)' )name,nelmts

 if(name==' 2nd derivatives (non-stat.)  - ' .or. name==' 2rd derivatives (non-stat.)  - ')then
   ddb%typ(iblok)=1
 else if(name==' 2nd derivatives (stationary) - ' .or. name==' 2rd derivatives (stationary) - ')then
   ddb%typ(iblok)=2
 else if(name==' 3rd derivatives              - ')then
   ddb%typ(iblok)=3
 else if(name==' Total energy                 - ')then
   ddb%typ(iblok)=0
 else if(name==' 1st derivatives              - ')then
   ddb%typ(iblok)=4
 else if(name==' 2nd eigenvalue derivatives   - ' .or. name==' 2rd eigenvalue derivatives   - ')then
   ddb%typ(iblok)=5
 else
   write(message,'(6a)')&
   'The following string appears in the DDB in place of',&
   ' the block type description :',ch10,trim(name),ch10,&
   'Action: check your DDB.'
   MSG_ERROR(message)
 end if

!Read the 2nd derivative block
 if(ddb%typ(iblok)==1.or.ddb%typ(iblok)==2)then

!  First check if there is enough space to read it
   if(msize<(3*mpert*3*mpert))then
     write(message,'(3a)')&
     'There is not enough space to read a second-derivative block.',ch10,&
     'Action: increase msize and recompile.'
     MSG_ERROR(message)
   end if

!  Read the phonon wavevector
   read(nunit, '(4x,3es16.8,f6.1)' )(ddb%qpt(ii,iblok),ii=1,3),ddb%nrm(1,iblok)

!  Read every element
   do ii=1,nelmts
     read(nunit,*)idir1,ipert1,idir2,ipert2,ar,ai
     index=idir1+3*((ipert1-1)+mpert*((idir2-1)+3*(ipert2-1)))
     ddb%flg(index,iblok)=1
     ddb%val(1,index,iblok)=ar
     ddb%val(2,index,iblok)=ai
   end do

!  Read the 3rd derivative block
 else if(ddb%typ(iblok)==3)then

!  First check if there is enough space to read it
   if(msize<(3*mpert*3*mpert*3*mpert))then
     write(message, '(a,a,a,i10,a,i10,a,a,a)' )&
     'There is not enough space to read a third-derivative block.',ch10,&
     'The size provided is only ',msize,' although ',3*mpert*3*mpert*3*mpert,' is needed.',ch10,&
     'Action: increase msize and recompile.'
     MSG_ERROR(message)
   end if

!  Read the perturbation wavevectors
   read(nunit,'(4x,3es16.8,f6.1)')(ddb%qpt(ii,iblok),ii=1,3),ddb%nrm(1,iblok)
   read(nunit,'(4x,3es16.8,f6.1)')(ddb%qpt(ii,iblok),ii=4,6),ddb%nrm(2,iblok)
   read(nunit,'(4x,3es16.8,f6.1)')(ddb%qpt(ii,iblok),ii=7,9),ddb%nrm(3,iblok)

!  Read every element
   do ii=1,nelmts
     read(nunit,'(6i4,2d22.14)')idir1,ipert1,idir2,ipert2,idir3,ipert3,ar,ai
     index=idir1+                     &
&     3*((ipert1-1)+mpert*((idir2-1)+ &
&     3*((ipert2-1)+mpert*((idir3-1)+3*(ipert3-1)))))
     ddb%flg(index,iblok)=1
     ddb%val(1,index,iblok)=ar
     ddb%val(2,index,iblok)=ai
   end do

!  Read the total energy
 else if(ddb%typ(iblok)==0)then

!  First check if there is enough space to read it
   if(msize<1)then
     write(message, '(3a,i0,3a)' )&
&     'There is not enough space to read a total energy block.',ch10,&
&     'The size provided is only ',msize,' although 1 is needed.',ch10,&
&     'Action: increase msize and recompile.'
     MSG_ERROR(message)
   end if

!  Read the total energy
   read(nunit,'(2d22.14)')ar,ai
   ddb%flg(1,iblok)=1
   ddb%val(1,1,iblok)=ar
   ddb%val(2,1,iblok)=ai

!  Read the 1st derivative block
 else if (ddb%typ(iblok) == 4) then

!  First check if there is enough space to read it
   if (msize < (3*mpert)) then
     write(message, '(3a,i0,a,i0,3a)' )&
     'There is not enough space to read a first-derivative block.',ch10,&
     'The size provided is only ',msize,' although ',3*mpert,' is needed.',ch10,&
     'Action: increase msize and recompile.'
     MSG_ERROR(message)
   end if

!  Read every element
   do ii=1,nelmts
     read(nunit,'(2i4,2d22.14)')idir1,ipert1,ar,ai
     index=idir1 + 3*(ipert1 - 1)
     ddb%flg(index,iblok)=1
     ddb%val(1,index,iblok)=ar
     ddb%val(2,index,iblok)=ai
   end do

!  Read the 2nd eigenvalue derivative block
 else if(ddb%typ(iblok)==5)then

!  First check if there is enough space to read it
   if(msize<(3*mpert*3*mpert))then
     write(message, '(3a,i0,a,i0,3a)' )&
     'There is not enough space to read a second-derivative block.',ch10,&
     'The size provided is only ',msize,' although ',3*mpert*3*mpert*mband*nkpt,' is needed.',ch10,&
     'Action: increase msize and recompile.'
     MSG_ERROR(message)
   end if

!  Read the phonon wavevector
   read(nunit, '(4x,3es16.8,f6.1)' )(ddb%qpt(ii,iblok),ii=1,3),ddb%nrm(1,iblok)

!  Read the K point and band
   if(present(blkval2).and.present(kpt))then
     do ikpt=1,nkpt
       read(nunit, '(9x,3es16.8)')(kpt(ii,ikpt),ii=1,3)
       do iband=1,mband
         read(nunit, '(6x,i3)') band
!        Read every element
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

end subroutine ddb_read_block
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/rdddb9
!! NAME
!! rdddb9
!!
!! FUNCTION
!! This routine reads the derivative database entirely,
!! for use in ppddb9, and performs some checks and symmetrisation
!! At the end, the whole DDB is in central memory, contained in the
!! array ddb%val(2,msize,ddb%nblok).
!! The information on it is contained in the four arrays
!!   ddb%flg(msize,ddb%nblok) : blok flag for each element
!!   ddb%qpt(9,ddb%nblok)  : blok wavevector (unnormalized)
!!   ddb%nrm(3,ddb%nblok)  : blok wavevector normalization
!!   ddb%typ(ddb%nblok)    : blok type
!!
!! INPUTS
!! atifc(natom) = atifc(ia) equals 1 if the analysis of ifc has to be done for atom ia; otherwise 0
!! ddbun = unit number for DDB io
!! dimekb=dimension of ekb (for the time being, only for norm- conserving psps)
!! iout=unit number for output of formatted data
!! filnam=name of input file
!! lmnmax=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!       =if useylm=0, max number of (l,n)   comp. over all type of psps
!! mband=maximum number of bands
!! mpert =maximum number of ipert
!! msize=maximum size of data blocks
!! msym =maximum number of symmetry elements in space group
!! natifc = number of atoms for which the analysis of ifc is done
!! natom = number of atoms
!! ntypat=number of atom types
!! usepaw= 0 for non paw calculation; =1 for paw calculation
!!
!! OUTPUT
!! acell(3)=length scales of cell (bohr)
!! amu(ntypat)=mass of the atoms (atomic mass unit)
!! ddb : ddb blok datatype
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
!! PARENTS
!!      m_ddb
!!
!! CHILDREN
!!      ddb_hdr_free,ddb_hdr_open_read
!!
!! SOURCE

subroutine rdddb9(acell,atifc,amu,ddb,ddbun,filnam,gmet,gprim,indsym,iout,&
& mband,mpert,msize,msym,natifc,natom,nkpt,nsym,ntypat,&
& rmet,rprim,symrec,symrel,symafm,tnons,typat,ucvol,xcart,xred,zion,znucl)

!Arguments ------------------------------------
! NOTE: these are used for dimensioning and then re-assigned in ioddb8.
!   This is almost definitely bad practice. In particular
!    it should be indsym(4,msym,natom),
!   and
!    the allocation allocate(kpt(3,nkpt)) is strange
!scalars
 integer,intent(in) :: ddbun,iout,mband,mpert,msize,msym,natifc
 integer,intent(inout) :: natom,nkpt,nsym,ntypat
 real(dp),intent(out) :: ucvol
 character(len=*),intent(in) :: filnam
 type(ddb_type),intent(inout) :: ddb
!arrays
 integer,intent(inout) :: atifc(natom)
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
 integer :: iblok,isym
 integer :: nsize,timrev
 integer :: i1dir,i1pert,i2dir,i2pert,i3dir,i3pert,npert
 real(dp),parameter :: tolsym8=tol8
 character(len=500) :: message
 type(ddb_hdr_type) :: ddb_hdr
!arrays
 integer :: symq(4,2,msym)
 integer,allocatable :: car3flg(:,:,:,:,:,:),carflg(:,:,:,:)
 integer,allocatable :: tmpflg(:,:,:,:,:,:),rfpert(:,:,:,:,:,:)
 real(dp) :: gprimd(3,3),qpt(3),rprimd(3,3)
 real(dp),allocatable :: d2cart(:,:,:,:,:),d3cart(:,:,:,:,:,:,:)
 real(dp),allocatable :: tmpval(:,:,:,:,:,:,:)

! *********************************************************************

 DBG_ENTER("COLL")

!Open the input derivative database file and read the header
 call ddb_hdr_open_read(ddb_hdr, filnam, ddbun, DDB_VERSION, msym=msym, mband=mband)

 !nkpt = ddb_hdr%nkpt
 !ntypat = ddb_hdr%ntypat
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

 call ddb_hdr_free(ddb_hdr)

!Compute different matrices in real and reciprocal space, also
!checks whether ucvol is positive.
 call mkrdim(acell,rprim,rprimd)
 call metric(gmet,gprimd,iout,rmet,rprimd,ucvol)

!Obtain reciprocal space primitive transl g from inverse trans of r
!(Unlike in abinit, gprim is used throughout ifc; should be changed, later)
 call matr3inv(rprim,gprim)

!Generate atom positions in cartesian coordinates
 call xred2xcart(natom,rprimd,xcart,xred)

!Transposed inversion of the symmetry matrices, for use in the reciprocal space
 do isym=1,nsym
   call mati3inv(symrel(:,:,isym),symrec(:,:,isym))
 end do

!SYMATM generates for all the atoms and all the symmetries, the atom
!on which the referenced one is sent and also the translation bringing
!back this atom to the referenced unit cell
 call symatm(indsym,natom,nsym,symrec,tnons,tolsym8,typat,xred)

!Check the correctness of some input parameters, and perform small treatment if needed.
 call chkin9(atifc,natifc,natom)

!Read the blocks from the input database, and close it.
 write(message, '(3a,i0,a)' )ch10,ch10,' rdddb9: read ',ddb%nblok,' blocks from the input DDB '
 call wrtout(std_out,message,'COLL')

 do iblok=1,ddb%nblok
   call ddb%read_block(iblok,mband,mpert,msize,nkpt,ddbun)

   !  Here complete the matrix by symmetrisation of the existing elements
   if(ddb%typ(iblok)==1 .or. ddb%typ(iblok)==2) then

     qpt(1)=ddb%qpt(1,iblok)/ddb%nrm(1,iblok)
     qpt(2)=ddb%qpt(2,iblok)/ddb%nrm(1,iblok)
     qpt(3)=ddb%qpt(3,iblok)/ddb%nrm(1,iblok)

     ! Examine the symmetries of the q wavevector
     call littlegroup_q(nsym,qpt,symq,symrec,symafm,timrev,prtvol=0)

     nsize=3*mpert*3*mpert
     ABI_MALLOC(tmpflg,(3,mpert,3,mpert,1,1))
     ABI_MALLOC(tmpval,(2,3,mpert,3,mpert,1,1))

     tmpflg(:,:,:,:,1,1) = reshape(ddb%flg(1:nsize,iblok), shape = (/3,mpert,3,mpert/))
     tmpval(1,:,:,:,:,1,1) = reshape(ddb%val(1,1:nsize,iblok), shape = (/3,mpert,3,mpert/))
     tmpval(2,:,:,:,:,1,1) = reshape(ddb%val(2,1:nsize,iblok), shape = (/3,mpert,3,mpert/))

     ! Then apply symmetry operations
     call d2sym3(tmpflg,tmpval,indsym,mpert,natom,nsym,qpt,symq,symrec,symrel,timrev,1)

     ! Transform the dynamical matrix in cartesian coordinates
     ABI_MALLOC(carflg,(3,mpert,3,mpert))
     ABI_MALLOC(d2cart,(2,3,mpert,3,mpert))

     call cart29(tmpflg,tmpval,carflg,d2cart,gprimd,1,mpert,natom,1,ntypat,rprimd,typat,ucvol,zion)

     ddb%flg(1:nsize,iblok) = reshape(carflg,shape = (/3*mpert*3*mpert/))
     ddb%val(1,1:nsize,iblok) = reshape(d2cart(1,:,:,:,:), shape = (/3*mpert*3*mpert/))
     ddb%val(2,1:nsize,iblok) = reshape(d2cart(2,:,:,:,:), shape = (/3*mpert*3*mpert/))

     ABI_FREE(carflg)
     ABI_FREE(d2cart)
     ABI_FREE(tmpflg)
     ABI_FREE(tmpval)

   else if (ddb%typ(iblok) == 3) then

     npert = 8
     nsize=3*npert*3*npert*3*npert
     ABI_MALLOC(tmpflg,(3,npert,3,npert,3,npert))
     ABI_MALLOC(tmpval,(2,3,npert,3,npert,3,npert))
     ABI_MALLOC(rfpert,(3,npert,3,npert,3,npert))

     tmpflg(:,:,:,:,:,:) = reshape(ddb%flg(1:nsize,iblok), shape = (/3,npert,3,npert,3,npert/))
     tmpval(1,:,:,:,:,:,:) = reshape(ddb%val(1,1:nsize,iblok), shape = (/3,npert,3,npert,3,npert/))
     tmpval(2,:,:,:,:,:,:) = reshape(ddb%val(2,1:nsize,iblok), shape = (/3,npert,3,npert,3,npert/))

!    Set the elements that are zero by symmetry for raman and
!    non-linear optical susceptibility tensors
     rfpert = 0
     rfpert(:,natom+2,:,natom+2,:,natom+2) = 1
     rfpert(:,1:natom,:,natom+2,:,natom+2) = 1
     rfpert(:,natom+2,:,1:natom,:,natom+2) = 1
     rfpert(:,natom+2,:,natom+2,:,1:natom) = 1
     call sytens(indsym,npert,natom,nsym,rfpert,symrec,symrel)
     do i1pert = 1,npert
       do i2pert = 1,npert
         do i3pert = 1,npert
           do i1dir=1,3
             do i2dir=1,3
               do i3dir=1,3
                 if ((rfpert(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)==-2) .and. &
&                  (tmpflg(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)/=1)) then
                   tmpval(:,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = zero
                   tmpflg(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)=1
                 end if
               end do
             end do
           end do
         end do
       end do
     end do

     call d3sym(tmpflg,tmpval,indsym,npert,natom,nsym,symrec,symrel)

     ABI_MALLOC(d3cart,(2,3,npert,3,npert,3,npert))
     ABI_MALLOC(car3flg,(3,npert,3,npert,3,npert))

     call nlopt(tmpflg,car3flg,tmpval,d3cart,gprimd,npert,natom,rprimd,ucvol)

     ddb%flg(1:nsize,iblok) = reshape(car3flg, shape = (/3*npert*3*npert*3*npert/))
     ddb%val(1,1:nsize,iblok) = reshape(d3cart(1,:,:,:,:,:,:), shape = (/3*npert*3*npert*3*npert/))
     ddb%val(2,1:nsize,iblok) = reshape(d3cart(2,:,:,:,:,:,:), shape = (/3*npert*3*npert*3*npert/))

     ABI_FREE(d3cart)
     ABI_FREE(car3flg)
     ABI_FREE(tmpflg)
     ABI_FREE(tmpval)
     ABI_FREE(rfpert)
   end if
 end do ! iblok

 close(ddbun)

 write(message,'(a)' )' Now the whole DDB is in central memory '
 call wrtout(std_out,message,'COLL')
 call wrtout(iout,message,'COLL')

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
!! PARENTS
!!      m_ddb,thmeig
!!
!! CHILDREN
!!      ddb_hdr_free,ddb_hdr_open_read
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
 character(len=500) :: message
!arrays
 integer,allocatable :: work(:)

! *********************************************************************

 if(natifc>natom)then
   write(message, '(a,i0,3a,i0,3a)' )&
&   'The number of atom ifc in the input files',natifc,',',ch10,&
&   'is larger than the number of atoms',natom,'.',ch10,&
&   'Action: change natifc in the input file.'
   MSG_ERROR(message)
 end if

 if(natifc>=1)then
   ABI_MALLOC(work,(natom))
   work(:)=0

   do iatifc=1,natifc
     if(atifc(iatifc)<=0.or.atifc(iatifc)>natom)then
       write(message, '(a,i0,5a,i0,3a)' )&
&       'For iatifc=',iatifc,', the number of the atom ifc to be ',ch10,&
&       'analysed is not valid : either negative, ',ch10,&
&       'zero, or larger than natom =',natom,'.',ch10,&
&       'Action: change atifc in your input file.'
       MSG_ERROR(message)
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
!! PARENTS
!!      m_ddb,nonlinear
!!
!! CHILDREN
!!      ddb_hdr_free,ddb_hdr_open_read
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
&             (blkflg(i1dir,i1pert,i3dir,i3pert,i2dir,i2pert)/=0).and. &
&             (blkflg(i2dir,i2pert,i1dir,i1pert,i3dir,i3pert)/=0).and. &
&             (blkflg(i2dir,i2pert,i3dir,i3pert,i1dir,i1pert)/=0).and. &
&             (blkflg(i3dir,i3pert,i1dir,i1pert,i2dir,i2pert)/=0).and. &
&             (blkflg(i3dir,i3pert,i2dir,i2pert,i1dir,i1pert)/=0)) then

               d3cart(:,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = &
&               (  d3(:,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) + &
&               d3(:,i1dir,i1pert,i3dir,i3pert,i2dir,i2pert) + &
&               d3(:,i2dir,i2pert,i1dir,i1pert,i3dir,i3pert) + &
&               d3(:,i2dir,i2pert,i3dir,i3pert,i1dir,i1pert) + &
&               d3(:,i3dir,i3pert,i1dir,i1pert,i2dir,i2pert) + &
&               d3(:,i3dir,i3pert,i2dir,i2pert,i1dir,i1pert))*sixth

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

!Compute non linear-optical coefficients d_ijk (atomic units)

 i1pert = natom+2
 d3cart(:,:,i1pert,:,i1pert,:,i1pert) = -3._dp*d3cart(:,:,i1pert,:,i1pert,:,i1pert)/(ucvol*2._dp)

!Compute first-order change in the electronic dielectric
!susceptibility (Bohr^-1) induced by an atomic displacement
 d3cart(1:2,1:3,1:natom,1:3,natom + 2,1:3,natom + 2) = &
& -6._dp*d3cart(1:2,1:3,1:natom,1:3,natom + 2,1:3,natom + 2)/ucvol

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
!!
!! INPUTS
!!  filename=DDB filename.
!!  brav = 1 or -1 -> simple lattice; 2 -> face-centered cubic;
!!         3 -> body-centered lattice; 4 -> hexagonal lattice (D6h)
!!  natom=Number of atoms in the unit cell
!!  atifc(natom)=list of the atom ifc to be analysed
!!  natifc = number of atoms for which the analysis of ifc is done
!!  comm=MPI communicator.
!!  [prtvol] = Verbosity level
!!
!! OUTPUT
!!  ddb<type(ddb_type)>=Object storing the DDB results.
!!  Crystal<type(crystal_t)>=Crystal structure parameters
!!  atifc(natom) =  atifc(ia) equals 1 if the analysis of ifc
!!    has to be done for atom ia; otherwise 0.
!!
!! TODO
!!   Sorry for the presence of natom, natifc and atifc.
!!   They are needed for legacy code!
!!
!! PARENTS
!!      anaddb,dfpt_looppert,eph,m_effective_potential_file,m_gruneisen
!!
!! CHILDREN
!!      ddb_hdr_free,ddb_hdr_open_read
!!
!! SOURCE

subroutine ddb_from_file(ddb, filename, brav, natom, natifc, atifc, crystal, comm, prtvol)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm,brav,natom,natifc
 integer,optional,intent(in) :: prtvol
 character(len=*),intent(in) :: filename
 type(crystal_t),intent(out) :: Crystal
 type(ddb_type),intent(inout) :: ddb
!array
 integer,intent(inout) :: atifc(natom)

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: ierr,ii,msym,dimekb,lmnmax,mband,nkpt,ntypat,nsym,usepaw
 integer :: mtyp,mpert,msize,ddb_natom,nblok,occopt,timrev,space_group,npsp,ddbun
 real(dp) :: factor,ucvol
 logical :: use_antiferro
 type(ddb_hdr_type) :: ddb_hdr
!arrays
 integer,allocatable :: symrec(:,:,:),symrel(:,:,:),symafm(:),indsym(:,:,:),typat(:)
 real(dp) :: acell(3),gmet(3,3),gprim(3,3),rmet(3,3),rprim(3,3),rprimd(3,3)
 real(dp),allocatable :: amu(:),xcart(:),xred(:,:),zion(:),znucl(:),tnons(:,:)
 character(len=132),allocatable :: title(:)
 character(len=500) :: message

! ************************************************************************

 DBG_ENTER("COLL")

! Must read natom from the DDB before being able to allocate some arrays needed for invars9
 ddbun = get_unit()
 call ddb_hdr_open_read(ddb_hdr, filename, ddbun, DDB_VERSION, comm=comm, dimonly=1)

 nblok = ddb_hdr%nblok
 mtyp = ddb_hdr%mblktyp
 msym = ddb_hdr%msym
 ddb_natom = ddb_hdr%natom
 ntypat = ddb_hdr%ntypat
 mband = ddb_hdr%mband
 nkpt = ddb_hdr%nkpt
 usepaw = ddb_hdr%usepaw
 dimekb = ddb_hdr%psps%dimekb
 lmnmax = ddb_hdr%psps%lmnmax

 ! JWZ occopt was used below before being initialized 13 April 2018
 occopt = ddb_hdr%occopt
 call ddb_hdr_free(ddb_hdr)

 if (ddb_natom /= natom) then
   MSG_ERROR(sjoin("input natom:",itoa(natom),"does not agree with DDB value:",itoa(natom)))
 end if

 mpert = natom+MPERT_MAX
 msize=3*mpert*3*mpert; if (mtyp==3) msize=msize*3*mpert

 ! Allocate arrays depending on msym (which is actually fixed to nsym inside inprep8)
 ABI_MALLOC(symrel,(3,3,msym))
 ABI_MALLOC(symafm,(msym))
 ABI_MALLOC(tnons,(3,msym))
 ABI_MALLOC(typat,(natom))
 ABI_MALLOC(xred,(3,natom))
 ABI_MALLOC(zion,(ntypat))
 ABI_MALLOC(znucl,(ntypat))

 ! Master reads and then broadcasts data.
 if (xmpi_comm_rank(comm) == master) then

   ABI_MALLOC(symrec,(3,3,msym))
   ABI_MALLOC(indsym,(4,msym,natom))
   ABI_MALLOC(xcart,(3*natom))
   ABI_MALLOC(amu,(ntypat))

   ddbun = get_unit() ! FIXME: The treatment of the unit number in rdddb9 is ugly!

   call ddb%malloc(msize, nblok, natom, ntypat)

   call rdddb9(acell,atifc,amu,ddb,&
    ddbun,filename,gmet,gprim,indsym,ab_out,&
    mband,mpert,msize,msym,&
    natifc,ddb_natom,nkpt,nsym,ntypat,&
    rmet,rprim,symrec,symrel,symafm,&
    tnons,typat,ucvol,xcart,xred,zion,znucl)

   close(ddbun)

   ABI_FREE(symrec)
   ABI_FREE(indsym)
   ABI_FREE(xcart)

   ! Renormalize rprim to possibly satisfy the constraint abs(rprim(1,2))=half when abs(brav)/=1
   ! This section is needed to preserver the behaviour of the old implementation.
   if (abs(brav)/=1 .and. abs(abs(rprim(1,2))-half)>tol10) then
     if(abs(rprim(1,2))<tol6)then
       write(message, '(a,i0,7a)' )&
        'The input DDB value of brav is ',brav,',',ch10,&
        'and the one of rprim(1,2) is zero.',ch10,&
        'These are incompatible',ch10,&
        'Action: check the value of brav and rprim(1,2) in your DDB.'
       MSG_ERROR(message)
     end if
     factor=abs(rprim(1,2))*two
     acell(:)=acell(:)*factor
     rprim(:,:)=rprim(:,:)/factor
     gprim(:,:)=gprim(:,:)*factor
   end if

   ! Save variables needed to call legacy code.
   ddb%acell = acell
   ddb%rprim = rprim
   ddb%gprim = gprim

   ! Other useful quantities.
   ! 2 is to preserve the old behaviour
   ddb%prtvol = 2; if (present(prtvol)) ddb%prtvol = prtvol
   ddb%occopt = occopt
   ddb%amu = amu
   ABI_FREE(amu)

   ! Now the whole DDB is in central memory, contained in the array ddb%val(2,msize,nblok).
   ! The data is contained in the four arrays
   !   ddb%flg(msize,nblok) : blok flag for each element
   !   ddb%qpt(9,nblok)     : blok wavevector (unnormalized)
   !   ddb%nrm(3,nblok)     : blok wavevector normalization
   !   ddb%typ(nblok)       : blok type
 end if

 if (xmpi_comm_size(comm) > 1) then
   call ddb%bcast(master, comm)
   call xmpi_bcast(atifc, master, comm, ierr)
   call xmpi_bcast(nsym, master, comm, ierr)
   call xmpi_bcast(symrel, master, comm, ierr)
   call xmpi_bcast(symafm, master, comm, ierr)
   call xmpi_bcast(typat, master, comm, ierr)
   call xmpi_bcast(acell, master, comm, ierr)
   call xmpi_bcast(occopt, master, comm, ierr)
   call xmpi_bcast(gprim, master, comm, ierr)
   call xmpi_bcast(rprim, master, comm, ierr)
   call xmpi_bcast(tnons, master, comm, ierr)
   call xmpi_bcast(xred, master, comm, ierr)
   call xmpi_bcast(zion, master, comm, ierr)
   call xmpi_bcast(znucl, master, comm, ierr)
 end if

 ! Initialize crystal_t object.
 call mkrdim(acell,rprim,rprimd)

!FIXME: These variables are hardcoded
 npsp = ntypat; space_group = 0; timrev = 2
 use_antiferro=.FALSE. !;  use_antiferro=(nspden==2.and.nsppol==1)
 ABI_MALLOC(title, (ntypat))

 do ii=1,ntypat
   write(title(ii),'(a,i0)')"No title for typat ",ii
 end do

 ! Warning znucl is dimensioned with ntypat = nspsp hence alchemy is not supported here
 call crystal_init(ddb%amu,Crystal,space_group,natom,npsp,ntypat,nsym,rprimd,typat,xred,&
   zion,znucl,timrev,use_antiferro,.FALSE.,title,&
   symrel=symrel(:,:,1:nsym),tnons=tnons(:,1:nsym),symafm=symafm(1:nsym))

 ABI_FREE(title)
 ABI_FREE(symrel)
 ABI_FREE(symafm)
 ABI_FREE(tnons)
 ABI_FREE(typat)
 ABI_FREE(xred)
 ABI_FREE(zion)
 ABI_FREE(znucl)

 DBG_EXIT("COLL")

end subroutine ddb_from_file
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
!! PARENTS
!!      thmeig
!!
!! CHILDREN
!!      ddb_hdr_free,ddb_hdr_open_read
!!
!! SOURCE

subroutine carttransf(blkflg,blkval2,carflg,gprimd,iqpt,mband,&
& mpert,msize,natom,nblok,nkpt,rprimd)

!Arguments ------------------------------------
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
integer :: iatom1,iatom2,iband,idir1,idir2,ikpt
integer :: index
!arrays
real(dp),allocatable :: blkflgtmp(:,:,:,:,:)
real(dp),allocatable :: blkval2tmp(:,:,:,:,:,:)
real(dp),allocatable :: d2cart(:,:,:,:,:)

! *********************************************************************

!Start by allocating local arrays
 ABI_MALLOC(blkflgtmp,(3,mpert,3,mpert,1))
 ABI_MALLOC(blkval2tmp,(2,3,mpert,3,mpert,1))
 ABI_MALLOC(d2cart,(2,3,mpert,3,mpert))

!Begin by formating the arrays to be compatible with cart29
!Then call cart29 to transform the arrays in cartesian coordinates
!Finally reformat the cartesian arrays in old format
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

!    The 1sin the argument of cart29 are respectively iblok and nblok. We are doing only one blok.
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

!Deallocating local arrays
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
!! PARENTS
!!      m_ddb
!!
!! CHILDREN
!!      ddb_hdr_free,ddb_hdr_open_read
!!
!! SOURCE

subroutine carteig2d(blkflg,blkval,carflg,d2cart,&
& gprimd,iblok,mpert,natom,nblok,rprimd)

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

!First, copy the data blok in place.
 d2cart(:,:,:,:,:)=blkval(:,:,:,:,:,iblok)

!Cartesian coordinates transformation (in two steps)
!First step
 do ipert1=1,mpert
   do ipert2=1,mpert
     do ii=1,2
       do idir1=1,3
         do idir2=1,3
           vec1(idir2)=d2cart(ii,idir1,ipert1,idir2,ipert2)
!          Note here blkflg
           flg1(idir2)=blkflg(idir1,ipert1,idir2,ipert2,iblok)
         end do
         call cart39(flg1,flg2,gprimd,ipert2,natom,rprimd,vec1,vec2)
         do idir2=1,3
           d2cart(ii,idir1,ipert1,idir2,ipert2)=vec2(idir2)
!          And here carflg
           carflg(idir1,ipert1,idir2,ipert2)=flg2(idir2)
         end do
       end do
     end do
   end do
 end do

!Second step
 do ipert1=1,mpert
   do ipert2=1,mpert
     do ii=1,2
       do idir2=1,3
         do idir1=1,3
           vec1(idir1)=d2cart(ii,idir1,ipert1,idir2,ipert2)
!          Note here carflg
           flg1(idir1)=carflg(idir1,ipert1,idir2,ipert2)
         end do
         call cart39(flg1,flg2,gprimd,ipert1,natom,rprimd,vec1,vec2)
         do idir1=1,3
           d2cart(ii,idir1,ipert1,idir2,ipert2)=vec2(idir1)
!          And here carflg again
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
!! PARENTS
!!      m_ddb
!!
!! CHILDREN
!!      ddb_hdr_free,ddb_hdr_open_read
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
!! PARENTS
!!      m_ddb
!!
!! CHILDREN
!!      ddb_hdr_free,ddb_hdr_open_read
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

!Extraction of non-linear optical coefficients
 do elfd1 = 1,3
   do elfd2 = 1,3
     do elfd3 = 1,3
       dchide(elfd1,elfd2,elfd3) = d3cart(1,elfd1,natom+2,elfd2,natom+2,elfd3,natom+2)
     end do
   end do
 end do

!Transform to Voigt notations
 voigtindex(:,1) = (/1,2,3,2,1,1/)
 voigtindex(:,2) = (/1,2,3,3,3,2/)
 do ivoigt = 1, 6
   elfd2 = voigtindex(ivoigt,1)
   elfd3 = voigtindex(ivoigt,2)
   do elfd1 = 1, 3
     dvoigt(elfd1,ivoigt) = 0.5_dp*(dchide(elfd1,elfd2,elfd3) + dchide(elfd1,elfd3,elfd2))
   end do
 end do

!Transform to pm/V
 dvoigt(:,:) = dvoigt(:,:)*16*(pi**2)*(Bohr_Ang**2)*1.0d-8*eps0/e_Cb

!Extraction of $\frac{d \chi}{d \tau}$
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
           sumrule(depl,elfd1,elfd2) = sumrule(depl,elfd1,elfd2) + &
&           dchidt(iatom,depl,elfd1,elfd2)
         end do
         do iatom = 1, natom
           dchidt(iatom,depl,elfd1,elfd2) = dchidt(iatom,depl,elfd1,elfd2) - &
&           wghtat(iatom)*sumrule(depl,elfd1,elfd2)
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
!!  iblock=Index of the block in the DDB file. 0 if not found.
!!
!! PARENTS
!!
!! CHILDREN
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
 rftyp = 0

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
!! PARENTS
!!
!! CHILDREN
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
!! Reads the Dielectric Tensor from the DDB file
!!
!! INPUTS
!!  ddb<type(ddb_type)>=Derivative database.
!!  rftyp  = 1 if non-stationary block
!!           2 if stationary block
!!           3 if third order derivatives
!!
!! OUTPUT
!!  dielt(3,3) = Macroscopic dielectric tensor
!!  iblok=Index of the block containing the data. 0 if block is not found.
!!
!! NOTES
!!  dielt is initialized to one_3D if the derivatives are not available in the DDB file.
!!
!! PARENTS
!!
!! CHILDREN
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
 character(len=1000) :: message
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

   write(message,'(a,3es16.6,3es16.6,3es16.6)' )&
&   ' Dielectric Tensor ',&
&   dielt(1,1),dielt(1,2),dielt(1,3),&
&   dielt(2,1),dielt(2,2),dielt(2,3),&
&   dielt(3,1),dielt(3,2),dielt(3,3)

   call wrtout(std_out,message,'COLL')

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
!! Reads the Dielectric Tensor from the DDB file
!!
!! INPUTS
!!  ddb<type(ddb_type)>=Derivative database.
!!  rftyp  = 1 if non-stationary block
!!           2 if stationary block
!!           3 if third order derivatives
!!
!! OUTPUT
!!  quadrupoles(3,3) = Macroscopic dielectric tensor
!!  iblok=Index of the block containing the data. 0 if block is not found.
!!
!! NOTES
!!  quadrupoles is initialized to zero if the derivatives are not available in the DDB file.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function ddb_get_quadrupoles(ddb, crystal, rftyp, quadrupoles) result(iblok)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: rftyp
 class(ddb_type),intent(in) :: ddb
 type(crystal_t),intent(in) :: crystal
!arrays
 real(dp),intent(out) :: quadrupoles(3,3,3,ddb%natom)

!Local variables -------------------------
!scalars
 integer :: ii, jj, iatdir, iatom, iq1dir, iq2dir
 integer :: quad_unt
 character(len=500) :: msg
!arrays
 integer :: rfelfd(4),rfphon(4),rfstrs(4)
 real(dp) :: qphnrm(3),qphon(3,3)

! *********************************************************************

 ! Look for the Gamma Block in the DDB
 qphon(:,:)=zero
 qphnrm(:)=one
 rfphon(1:2)=1
 rfelfd(1:2)=1
 rfstrs(1:2)=8

 call ddb%get_block(iblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,rftyp)

 ! Read the dielectric tensor only if the Gamma-block was found in the DDB
 ! In case it was not found, iblok = 0
 iblok = 0
 quadrupoles=zero

 !Temporary hack to read the quadrupole tensor from a text file
 if (.not.file_exists("quadrupoles_cart.out")) return
 quad_unt = 71
 open(unit=quad_unt,file="quadrupoles_cart.out",action="read")
 do ii=1,2
   read(quad_unt,*) msg
 end do

 do ii=1,3
   do jj=1,3*3*crystal%natom
     read(quad_unt,'(4(i5,3x),2(1x,f20.10))') iq2dir,iatom,iatdir,iq1dir,quadrupoles(iq1dir,iq2dir,iatdir,iatom)
     !write(*,*) iq2dir,iatom,iatdir,iq1dir,quadrupoles(iq1dir,iq2dir,iatdir,iatom)
   end do
   read(quad_unt,'(a)') msg
 end do
 close(quad_unt)
 iblok = 1

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
!! PARENTS
!!
!! CHILDREN
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
 rftyp = 3

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
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

type(asrq0_t) function ddb_get_asrq0(ddb, asr, rftyp, xcart) result(asrq0)

!Arguments ------------------------------------
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
   MSG_ERROR(sjoin("Wrong value for asr:", itoa(asr)))
 end select

end function ddb_get_asrq0
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
!! PARENTS
!!
!! CHILDREN
!!      ddb_hdr_free,ddb_hdr_open_read
!!
!! SOURCE

subroutine ddb_diagoq(ddb, crystal, qpt, asrq0, symdynmat, rftyp, phfrq, displ_cart, &
                      out_eigvec,out_displ_red)   ! Optional [out]

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: rftyp,symdynmat
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
 rfphon(1:2)=1
 rfelfd(1:2)=0
 rfstrs(1:2)=0
 qphon_padded = zero
 qphon_padded(:,1) = qpt
 natom = crystal%natom

 call ddb%get_block(iblok,qphon_padded,qphnrm,rfphon,rfelfd,rfstrs,rftyp)
 if (iblok == 0) then
   MSG_ERROR(sjoin("Cannot find q-point ", ktoa(qpt)," in DDB file"))
 end if

 ! Copy the dynamical matrix in d2cart
 d2cart(:,1:ddb%msize) = ddb%val(:,:,iblok)

 ! Eventually impose the acoustic sum rule based on previously calculated d2asr
 call asrq0%apply(natom, ddb%mpert, ddb%msize, crystal%xcart, d2cart)

 ! Calculation of the eigenvectors and eigenvalues of the dynamical matrix
 call dfpt_phfrq(ddb%amu,displ_cart,d2cart,eigval,eigvec,crystal%indsym,&
&  ddb%mpert,crystal%nsym,natom,crystal%nsym,crystal%ntypat,phfrq,qphnrm(1),my_qpt,&
&  crystal%rprimd,symdynmat,crystal%symrel,crystal%symafm,crystal%typat,crystal%ucvol)

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
!! PARENTS
!!      anaddb,m_ddb,m_phonons
!!
!! CHILDREN
!!      ddb_hdr_free,ddb_hdr_open_read
!!
!! SOURCE

subroutine asrq0_apply(asrq0, natom, mpert, msize, xcart, d2cart)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom, msize, mpert
 class(asrq0_t),intent(inout) :: asrq0
!arrays
 real(dp),intent(in) :: xcart(3,natom)
 real(dp),intent(inout) :: d2cart(2,msize)

! ************************************************************************

 if (asrq0%asr /= 0 .and. asrq0%iblok == 0) then
   MSG_WARNING("asr != 0 but DDB file does not contain q=Gamma. D(q) cannot be corrected")
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
   MSG_ERROR(sjoin("Wrong value for asr:", itoa(asrq0%asr)))
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
!! PARENTS
!!      anaddb,m_effective_potential_file
!!
!! CHILDREN
!!      ddb_hdr_free,ddb_hdr_open_read
!!
!! SOURCE

subroutine asrq0_free(asrq0)

!Arguments ------------------------------------
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

!!****f* m_ddb/ddb_write_block
!!
!! NAME
!! ddb_write_block
!!
!! FUNCTION
!! This routine writes blocks of data in the DDBs.
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
!! PARENTS
!!      ddb_interpolate,dfptnl_doutput,gstate,mblktyp1,mblktyp5
!!
!! CHILDREN
!!      ddb_hdr_free,ddb_hdr_open_read
!!
!! SOURCE

subroutine ddb_write_block(ddb,iblok,choice,mband,mpert,msize,nkpt,nunit,&
&     blkval2,kpt) !optional

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

! *********************************************************************


!Count the number of elements
 nelmts=0
 do ii=1,msize
   if(ddb%flg(ii,iblok)==1)nelmts=nelmts+1
 end do

!Write the block type and number of elements
 write(nunit,*)' '
 if (ddb%typ(iblok) == 0) then
   write(nunit, '(a,i8)' )' Total energy                 - # elements :',nelmts
 else if (ddb%typ(iblok)==1) then
   write(nunit, '(a,i8)' )' 2nd derivatives (non-stat.)  - # elements :',nelmts
 else if(ddb%typ(iblok)==2) then
   write(nunit, '(a,i8)' )' 2nd derivatives (stationary) - # elements :',nelmts
 else if(ddb%typ(iblok)==3) then
   write(nunit, '(a,i8)' )' 3rd derivatives              - # elements :',nelmts
 else if (ddb%typ(iblok) == 4) then
   write(nunit, '(a,i8)' )' 1st derivatives              - # elements :',nelmts
 else if (ddb%typ(iblok) == 5) then
   write(nunit, '(a,i8)' )' 2nd eigenvalue derivatives   - # elements :',nelmts
 end if

!Write the 2nd derivative block
 if(ddb%typ(iblok)==1.or.ddb%typ(iblok)==2)then

!  Write the phonon wavevector
   write(nunit, '(a,3es16.8,f6.1)' )' qpt',(ddb%qpt(ii,iblok),ii=1,3),ddb%nrm(1,iblok)

!  Write the matrix elements
   if(choice==2)then
     ii=0
     do ipert2=1,mpert
       do idir2=1,3
         do ipert1=1,mpert
           do idir1=1,3
             ii=ii+1
             if(ddb%flg(ii,iblok)==1)then
               write(nunit,'(4i4,2d22.14)')idir1,ipert1,idir2,ipert2,ddb%val(1,ii,iblok),ddb%val(2,ii,iblok)
             end if
           end do
         end do
       end do
     end do
   end if

!  Write the 3rd derivative block
 else if(ddb%typ(iblok)==3)then

!  Write the phonon wavevectors
   write(nunit, '(a,3es16.8,f6.1)' )' qpt',(ddb%qpt(ii,iblok),ii=1,3),ddb%nrm(1,iblok)
   write(nunit, '(a,3es16.8,f6.1)' )'    ',(ddb%qpt(ii,iblok),ii=4,6),ddb%nrm(2,iblok)
   write(nunit, '(a,3es16.8,f6.1)' )'    ',(ddb%qpt(ii,iblok),ii=7,9),ddb%nrm(3,iblok)

!  Write the matrix elements
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
&                   idir1,ipert1,idir2,ipert2,idir3,ipert3,ddb%val(1,ii,iblok),ddb%val(2,ii,iblok)
                 end if
               end do
             end do
           end do
         end do
       end do
     end do
   end if

!  Write total energy
 else if (ddb%typ(iblok) == 0) then
   if (choice == 2) then
     write(nunit,'(2d22.14)')ddb%val(1,1,iblok),ddb%val(2,1,iblok)
   end if

!  Write the 1st derivative blok
 else if (ddb%typ(iblok) == 4) then
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

 else if (ddb%typ(iblok)==5) then
!  Write the phonon wavevector
   write(nunit, '(a,3es16.8,f6.1)' )' qpt',(ddb%qpt(ii,iblok),ii=1,3),ddb%nrm(1,iblok)
!  Write the matrix elements
   if(choice==2)then
     if(present(blkval2).and.present(kpt))then
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
                     write(nunit,'(4i4,2d22.14)')idir1,ipert1,idir2,ipert2,&
&                     blkval2(1,ii,iband,ikpt),blkval2(2,ii,iband,ikpt)
                   end if
                 end do !idir1
               end do  !ipert1
             end do   !idir2
           end do    !ipert2
         end do     !iband
       end do      !ikpt
     end if !blkval2
   end if !choice
 end if !ddb%typ(iblok)

end subroutine ddb_write_block
!!***

!!****f* m_ddb/dfptnl_doutput
!! NAME
!! dfptnl_doutput
!!
!! FUNCTION
!! Write the matrix of third-order derivatives to the output file and the DDB
!!
!! INPUTS
!!  blkflg(3,mpert,3,mpert,3,mpert)= ( 1 if the element of the 3dte
!!   has been calculated ; 0 otherwise )
!!  d3(2,3,mpert,3,mpert,3,mpert)= matrix of the 3DTE
!!  mpert =maximum number of ipert
!!  natom=Number of atoms
!!  ntypat=Number of type of atoms
!!  unddb = unit number for DDB output
!!
!! NOTES
!!  d3 holds the third-order derivatives before computing
!!  the permutations of the perturbations.
!!
!! PARENTS
!!      nonlinear
!!
!! CHILDREN
!!      ddb_free,ddb_malloc,ddb_write_block,wrtout
!!
!! SOURCE

subroutine dfptnl_doutput(blkflg,d3,mband,mpert,nkpt,natom,ntypat,unddb)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: mband,mpert,nkpt,unddb,natom,ntypat
!arrays
 integer,intent(in) :: blkflg(3,mpert,3,mpert,3,mpert)
 real(dp),intent(in) :: d3(2,3,mpert,3,mpert,3,mpert)

!Local variables -------------------------
!scalars
 integer :: choice,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert,index,msize
 character(len=500) :: message
 type(ddb_type) :: ddb

!*************************************************************************

 msize = 27*mpert*mpert*mpert
 call ddb%malloc(msize, 1, natom, ntypat)

 choice = 2

 ddb%typ = 3
 ddb%nrm = one
 ddb%qpt = zero   ! this has to be changed in case anharmonic
!force constants have been computed


!Write blok of third-order derivatives to ouput file

 write(message,'(a,a,a,a,a)')ch10,&
& ' Matrix of third-order derivatives (reduced coordinates)',ch10,&
& ' before computing the permutations of the perturbations',ch10
 call wrtout(ab_out,message,'COLL')

 write(ab_out,*)'    j1       j2       j3              matrix element'
 write(ab_out,*)' dir pert dir pert dir pert           real part           imaginary part'

 do i1pert=1,mpert
   do i1dir=1,3
     do i2pert=1,mpert
       do i2dir=1,3
         do i3pert=1,mpert
           do i3dir=1,3

             index = i1dir + &
&             3*((i1pert-1)+mpert*((i2dir-1) + &
&             3*((i2pert-1)+mpert*((i3dir-1) + 3*(i3pert-1)))))
             ddb%flg(index,1) = blkflg(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)
             ddb%val(:,index,1)= d3(:,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)

             if (blkflg(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)/=0) then

               write(ab_out,'(3(i4,i5),2f22.10)')&
&               i1dir,i1pert,i2dir,i2pert,i3dir,i3pert,&
&               d3(:,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)

             end if

           end do
         end do
       end do
     end do
   end do
 end do

!Write blok of third-order derivatives to DDB
 call ddb%write_block(1,choice,mband,mpert,msize,nkpt,unddb)

 call ddb_free(ddb)

end subroutine dfptnl_doutput
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_to_dtset
!! NAME
!! ddb_to_dtset
!!
!! FUNCTION
!!   Initialize a dataset object from ddb.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      ddb_hdr_free,ddb_hdr_open_read
!!
!! SOURCE


subroutine ddb_to_dtset(comm,dtset,filename,psps)

 !Arguments ------------------------------------
 integer,intent(in) :: comm
 type(dataset_type),intent(inout) :: dtset
 type(pseudopotential_type),intent(inout) :: psps
 ! type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)
 character(len=*),intent(in) :: filename
 !Local variables -------------------------
 integer :: mxnimage,ddbun
!integer :: ii, nn
 type(ddb_hdr_type) :: ddb_hdr

! ************************************************************************

 ABI_UNUSED(psps%usepaw)

!Set variables
 mxnimage = 1 ! Only 1 image in the DDB

! Must read natom from the DDB before being able to allocate some arrays needed for invars9
 ddbun = get_unit()
 call ddb_hdr_open_read(ddb_hdr,filename,ddbun,DDB_VERSION,comm=comm)
!close ddb file, just want to read the headers
 close(ddbun)
 dtset%ngfft = ddb_hdr%ngfft

! call psps_copy(psps, ddb_hdr%psps)

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
 if (allocated(dtset%acell_orig)) then
   ABI_DEALLOCATE(dtset%acell_orig)
 end if
 ABI_ALLOCATE(dtset%acell_orig,(3,mxnimage))
 dtset%acell_orig(1:3,1) = ddb_hdr%acell(:)

 if (allocated(dtset%rprim_orig)) then
   ABI_DEALLOCATE(dtset%rprim_orig)
 end if
 ABI_ALLOCATE(dtset%rprim_orig,(3,3,mxnimage))
 dtset%rprim_orig(1:3,1:3,1) = ddb_hdr%rprim(:,:)

 if (allocated(dtset%rprimd_orig)) then
   ABI_DEALLOCATE(dtset%rprimd_orig)
 end if
 ABI_ALLOCATE(dtset%rprimd_orig,(3,3,mxnimage))
 dtset%rprimd_orig(:,1,1) = ddb_hdr%rprim(:,1) * dtset%acell_orig(1,1)
 dtset%rprimd_orig(:,2,1) = ddb_hdr%rprim(:,2) * dtset%acell_orig(2,1)
 dtset%rprimd_orig(:,3,1) = ddb_hdr%rprim(:,3) * dtset%acell_orig(3,1)

 if (allocated(dtset%amu_orig)) then
   ABI_DEALLOCATE(dtset%amu_orig)
 end if
 ABI_ALLOCATE(dtset%amu_orig,(dtset%ntypat,mxnimage))
 dtset%amu_orig(:,1) = ddb_hdr%amu(:)

 if (allocated(dtset%typat)) then
   ABI_DEALLOCATE(dtset%typat)
 end if
 ABI_ALLOCATE(dtset%typat,(dtset%natom))
 dtset%typat(:) = ddb_hdr%typat(1:ddb_hdr%matom)

 if (allocated(dtset%spinat)) then
   ABI_DEALLOCATE(dtset%spinat)
 end if
 ABI_ALLOCATE(dtset%spinat,(3,dtset%natom))
 dtset%spinat(:,:) = ddb_hdr%spinat(1:3,1:ddb_hdr%matom)

 if (allocated(dtset%xred_orig)) then
   ABI_DEALLOCATE(dtset%xred_orig)
 end if
 ABI_ALLOCATE(dtset%xred_orig,(3,dtset%natom,mxnimage))
 dtset%xred_orig(:,:,1) = ddb_hdr%xred(1:3,1:ddb_hdr%matom)

 if (allocated(dtset%ziontypat)) then
   ABI_DEALLOCATE(dtset%ziontypat)
 end if
 ABI_ALLOCATE(dtset%ziontypat,(dtset%ntypat))
 dtset%ziontypat(1:ddb_hdr%mtypat) = ddb_hdr%zion(1:ddb_hdr%mtypat)

 if (allocated(dtset%znucl)) then
   ABI_DEALLOCATE(dtset%znucl)
 end if
 ABI_ALLOCATE(dtset%znucl,(dtset%ntypat))
 dtset%znucl(:) = ddb_hdr%znucl(1:ddb_hdr%mtypat)

 if (allocated(dtset%nband)) then
   ABI_DEALLOCATE(dtset%nband)
 end if
 ABI_ALLOCATE(dtset%nband,(dtset%nkpt))
 dtset%nband(:) = ddb_hdr%nband(1:ddb_hdr%mkpt*ddb_hdr%nsppol)

 if (allocated(dtset%symafm)) then
   ABI_DEALLOCATE(dtset%symafm)
 end if
 ABI_ALLOCATE(dtset%symafm,(dtset%nsym))
 dtset%symafm(:) = ddb_hdr%symafm(1:ddb_hdr%msym)

 if (allocated(dtset%symrel)) then
   ABI_DEALLOCATE(dtset%symrel)
 end if
 ABI_ALLOCATE(dtset%symrel,(3,3,dtset%nsym))
 dtset%symrel(:,:,:) = ddb_hdr%symrel(1:3,1:3,1:ddb_hdr%msym)

 if (allocated(dtset%tnons)) then
   ABI_DEALLOCATE(dtset%tnons)
 end if
 ABI_ALLOCATE(dtset%tnons,(3,dtset%nsym))
 dtset%tnons(:,:) = ddb_hdr%tnons(1:3,1:ddb_hdr%msym)

 if (allocated(dtset%kpt)) then
   ABI_DEALLOCATE(dtset%kpt)
 end if
 ABI_ALLOCATE(dtset%kpt,(3,dtset%nkpt))
 dtset%kpt(:,:) = ddb_hdr%kpt(1:3,1:ddb_hdr%mkpt)

 if (allocated(dtset%wtk)) then
   ABI_DEALLOCATE(dtset%wtk)
 end if
 ABI_ALLOCATE(dtset%wtk,(dtset%nkpt))
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

 call ddb_hdr_free(ddb_hdr)

end subroutine ddb_to_dtset
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/mblktyp1
!!
!! NAME
!! mblktyp1
!!
!! FUNCTION
!! This routine merges the derivative databases of type 0-4:
!! Total energy, (2nd derivatives (non-stat.),2nd derivatives (stationary),
!! 3rd derivatives, 1st derivatives
!!
!! The heading of the database is read, then the heading
!! of the temporary database to be added is read,
!! the code check their compatibility, and create a new
!! database that mixes the old and the temporary ones.
!! This process can be iterated.
!! The whole database will be stored in central memory.
!!
!! INPUTS
!!     chkopt=option for consistency checks between DDB files
!!     ddbun=define input and output unit numbers
!!     dscrpt=description of the output file
!!     filnam=name of input or output file
!!     mddb=maximum number of databases (cannot be made dynamic)
!!     nddb=number of input DDBs
!!     vrsddb=current version of the DDB
!!
!! OUTPUT
!!     msym=maximum number of symmetry elements in space group
!!     Merge the file
!!
!! PARENTS
!!      mrgddb
!!
!! CHILDREN
!!      ddb_free,ddb_hdr_compare,ddb_hdr_free,ddb_hdr_open_read
!!      ddb_hdr_open_write,ddb_malloc,ddb_write_block,read_block,wrtout
!!
!! SOURCE

subroutine mblktyp1(chkopt,ddbun,dscrpt,filnam,mddb,msym,nddb,vrsddb)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: chkopt,ddbun,mddb,nddb,vrsddb
 integer,intent(out) :: msym
 character(len=fnlen),intent(in) :: dscrpt
 character(len=fnlen),intent(in) :: filnam(mddb+1)

!Local variables -------------------------
!scalars
 integer :: choice,dimekb,iblok,iblok1,iblok2
 integer :: iddb,ii,lmnmax,matom
 integer :: mband,mblktyp,mblok,mkpt,mpert,msize,mtypat
 integer :: nblok,nblokt,nq
 integer :: tmerge,usepaw
 real(dp),parameter :: qtol=2.0d-8
 real(dp) :: diff
 type(ddb_type) :: ddb
 type(ddb_hdr_type) :: ddb_hdr, ddb_hdr8
 character(len=500) :: message
!arrays
 integer,allocatable :: mgblok(:)!,lloc(:)

! *********************************************************************

 ! Make sure there is more than one ddb to be read
 if(nddb==1)then
   write(message, '(a,a,a,a,a)' )&
&   'The initialisation mode of MRGDDB, that uses nddb=1,',&
&   'has been disabled in version 2.1 of ABINIT.',&
&   'Action: you should use DDBs that include the symmetry',&
&   'information (and that can be used and merged without',&
&   'initialisation), or you should use ABINITv2.0.'
   MSG_ERROR(message)
 end if

!Evaluate the maximal dimensions of arrays
 dimekb=0 ; matom=0 ; mband=0  ; mblok=0 ; mkpt=0
 msize=0  ; mtypat=0 ; lmnmax=0 ; usepaw=0 ; mblktyp = 1
 msym=192

 do iddb=1,nddb
   call ddb_hdr_open_read(ddb_hdr, filnam(iddb+1), ddbun, vrsddb, dimonly=1)

   mblok=mblok+ddb_hdr%nblok
   mblktyp=max(mblktyp,ddb_hdr%mblktyp)
   matom=max(matom,ddb_hdr%matom)
   mkpt=max(mkpt,ddb_hdr%mkpt)
   mtypat=max(mtypat,ddb_hdr%mtypat)
   msym=max(msym,ddb_hdr%msym)
   mband=max(mband,ddb_hdr%mband)
   dimekb=max(dimekb,ddb_hdr%psps%dimekb)
   lmnmax=max(lmnmax,ddb_hdr%psps%lmnmax)
   usepaw=max(usepaw,ddb_hdr%usepaw)

   call ddb_hdr_free(ddb_hdr)
 end do

 mpert=matom+MPERT_MAX
 msize=3*mpert*3*mpert
 if(mblktyp==3)msize=msize*3*mpert

 call ddb%malloc(msize, mblok, matom, mtypat)

!Allocate arrays
 ABI_ALLOCATE(mgblok,(mblok))

!**********************************************************************

!Read the first database

 write(std_out,*)' read the input derivative database information'
 call ddb_hdr_open_read(ddb_hdr, filnam(2), ddbun, vrsddb, &
& matom=matom,mtypat=mtypat,mband=mband,mkpt=mkpt,&
& msym=msym,dimekb=dimekb,lmnmax=lmnmax,usepaw=usepaw)

 if(ddb_hdr%nblok>=1)then
!  Read the blocks from the input database.
   write(message, '(a,i5,a)' ) ' read ',ddb_hdr%nblok, ' blocks from the input DDB '
   call wrtout(std_out,message,'COLL')
   do iblok=1,ddb_hdr%nblok
     call ddb%read_block(iblok,ddb_hdr%nband(1),mpert,msize,ddb_hdr%nkpt,ddbun)
!    Setup merged indicator
     mgblok(iblok)=0
   end do
 else
   write(message, '(a)' )' No bloks in the first ddb '
   call wrtout(std_out,message,'COLL')
 end if

!Close the first ddb
 close(ddbun)

!*********************************************

 nblok = ddb_hdr%nblok
!In case of merging of DDBs, iterate the reading
 do iddb=2,nddb

!  Open the corresponding input DDB, and read the database file information
   write(message, '(a,a,i6)' )ch10,' read the input derivative database number',iddb
   call wrtout(std_out,message,'COLL')

   call ddb_hdr_open_read(ddb_hdr8, filnam(iddb+1), ddbun, vrsddb, &
&   matom=matom,mtypat=mtypat,mband=mband,mkpt=mkpt,&
&   msym=msym,dimekb=dimekb,lmnmax=lmnmax,usepaw=usepaw)

   if (chkopt==1)then
!    Compare the current DDB and input DDB information.
!    In case of an inconsistency, halt the execution.
     call wrtout(std_out, ' compare the current and input DDB information', 'COLL')
     call ddb_hdr_compare(ddb_hdr, ddb_hdr8)

   else if(chkopt==0)then
!    No comparison between the current DDB and input DDB information.
     write(message, '(a)' )' no comparison between the current and input DDB information'
     call wrtout(std_out,message,'COLL')
     write(message, '(5a)' )&
&     'No comparison/check is performed for the current and input DDB information ',ch10,&
&     'because argument --nostrict was passed to the command line. ',ch10,&
&     'Use at your own risk!'
     MSG_COMMENT(message)
   end if

   call wrtout(std_out,' Will try to merge this input DDB with the current one.','COLL')

!  First estimate of the total number of bloks, and sto with error message if too large.
   write(message, '(a,i5)' ) ' Current number of bloks =',nblok
   call wrtout(std_out,message,'COLL')
   write(message, '(a,i5,a)' )' Will read ',ddb_hdr8%nblok,' blocks from the input DDB '
   call wrtout(std_out,message,'COLL')
   nblokt=nblok+ddb_hdr8%nblok
   if(nblokt>mblok)then
     write(message, '(a,i0,3a,i0,a)' )&
     'The expected number of blocks',nblokt,' is larger than',ch10,&
     'the maximum number of blocks',mblok,'.'
     MSG_ERROR(message)
   end if

!  Read the bloks from the temporary database, and close it.
!  Also setup the merging indicator
   do iblok=nblok+1,nblokt
     call ddb%read_block(iblok,ddb_hdr8%nband(1),mpert,msize,ddb_hdr8%nkpt,ddbun)
     mgblok(iblok)=0
   end do
   close(ddbun)

   nblok=nblokt
   write(message, '(a,i5)' ) ' Now, current number of bloks =',nblok
   call wrtout(std_out,message,'COLL')

   ! In certain cases, the different DDB will have different information
   ! on the pseudos (depending on fullinit)
   ! Here, we copy the information of the last DDB file,
   ! only to make the tests pass...
   ddb_hdr%psps%indlmn(:,:,:) = ddb_hdr8%psps%indlmn(:,:,:)
   ddb_hdr%psps%pspso(:) = ddb_hdr8%psps%pspso(:)
   ddb_hdr%psps%ekb(:,:) = ddb_hdr8%psps%ekb(:,:)

   call ddb_hdr_free(ddb_hdr8)
 end do

 call wrtout(std_out,' All DDBs have been read ','COLL')

!*********************************************************

!Check the equality of blocks, and eventually merge them

 if (nblok>=1) then
   call wrtout(std_out,' check the equality of blocks, and eventually merge ','COLL')
   do iblok2=2,nblok
     do iblok1=1,iblok2-1
       tmerge=0

!      Check the block type identity
       if(ddb%typ(iblok1)==ddb%typ(iblok2))then

!        Check the wavevector identities
         tmerge=1
         if(ddb%typ(iblok1)==1.or.ddb%typ(iblok1)==2)then
           nq=1
         else if(ddb%typ(iblok1)==3)then
!          Note : do not merge permutation related elements ....
           nq=3
         else if(ddb%typ(iblok1)==4 .or. ddb%typ(iblok1)==0)then
           nq=0
         end if
         if(nq/=0)then
           do ii=1,nq
             diff= ddb%qpt(1+3*(ii-1),iblok1)/ddb%nrm(ii,iblok1) - ddb%qpt(1+3*(ii-1),iblok2)/ddb%nrm(ii,iblok2)
             if (abs(diff) > qtol) tmerge=0
             diff=ddb%qpt(2+3*(ii-1),iblok1)/ddb%nrm(ii,iblok1) - ddb%qpt(2+3*(ii-1),iblok2)/ddb%nrm(ii,iblok2)
             if (abs(diff) > qtol) tmerge=0
             diff=ddb%qpt(3+3*(ii-1),iblok1)/ddb%nrm(ii,iblok1) - ddb%qpt(3+3*(ii-1),iblok2)/ddb%nrm(ii,iblok2)
             if (abs(diff) > qtol) tmerge=0
           end do ! ii
         end if

!        Now merges,
         if (tmerge == 1) then
           write(message, '(a,i5,a,i5)' )' merge block #',iblok2,' to block #',iblok1
           call wrtout(std_out,message,'COLL')
           mgblok(iblok2)=1
           do ii=1,msize
             if(ddb%flg(ii,iblok2)==1)then
               ddb%flg(ii,iblok1)=1
               ddb%val(1,ii,iblok1)=ddb%val(1,ii,iblok2)
               ddb%val(2,ii,iblok1)=ddb%val(2,ii,iblok2)
             end if
           end do
         end if

       end if
     end do
   end do

!  Count the final number of bloks
   tmerge=0
   do ii=1,nblok
     if(mgblok(ii)==1)tmerge=tmerge+1
   end do
   nblok=nblok-tmerge

!  Summarize the merging phase
   write(message, '(i6,a,i6,a)' )tmerge,' blocks are merged; the new DDB will have ',nblok,' blocks.'
   call wrtout(std_out,message,'COLL')
 end if !  End the condition on existence of more than one blok in current DDB

!**********************************************************************

 write(message, '(a,a)' )' open the output database, write the',' preliminary information '
 call wrtout(std_out,message,'COLL')

 ddb_hdr%dscrpt = trim(dscrpt)
 ddb_hdr%nblok = nblok
 ddb_hdr%mblktyp = mblktyp

 call ddb_hdr_open_write(ddb_hdr, filnam(1), ddbun, fullinit=1)

 if(nddb>1)then

!  Write the whole database
   call wrtout(std_out,' write the DDB ','COLL')
   choice=2
   do iblok=1,nblok+tmerge
     if(mgblok(iblok)==0)then
       write(std_out,'(a,i4)' ) ' Write bloc number',iblok
       call ddb%write_block(iblok,choice,ddb_hdr%nband(1),mpert,msize,ddb_hdr%nkpt,ddbun)
     else
       write(message, '(a,i4,a)' )&
&       ' Bloc number',iblok,' was merged, so do not write it'
       call wrtout(std_out,message,'COLL')
     end if
   end do

!  Also write summary of bloks at the end
   write(ddbun, '(/,a)' )' List of bloks and their characteristics '
   choice=3
   do iblok=1,nblok+tmerge
     if(mgblok(iblok)==0)then
       call ddb%write_block(iblok,choice,ddb_hdr%nband(1),mpert,msize,ddb_hdr%nkpt,ddbun)
     end if
   end do

 end if

 close (ddbun)

!*********************************************************************

!Deallocate arrays

 ABI_DEALLOCATE(mgblok)

 call ddb_hdr_free(ddb_hdr)
 call ddb_free(ddb)

end subroutine mblktyp1
!!***

!!****f* m_ddb/mblktyp5
!!
!! NAME
!! mblktyp5
!!
!! FUNCTION
!! This routine merges the derivative databases of type 5:
!! second-order eigenvalue derivatives
!!   why is this separate from mblktyp1? Should be merged at some point for consistency
!!
!! NOTES
!! The heading of the database is read, then the heading
!! of the temporary database to be added is read,
!! the code check their compatibility, and create a new
!! database that mixes the old and the temporary ones.
!!
!! TODO
!!  This routine is deprecated and will be removed
!!
!! INPUTS
!!     chkopt=option for consistency checks between DDB files
!!     codename=MRGDDB
!!     ddbun=define input and output unit numbers
!!     dscrpt=description of the output file
!!     filnam=name of input or output file
!!     mddb=maximum number of databases (cannot be made dynamic)
!!     nddb=number of input DDBs
!!
!! OUTPUT
!!     msym=maximum number of symmetry elements in space group
!!     Merge the file
!!
!! PARENTS
!!      mrgddb
!!
!! CHILDREN
!!      ddb_free,ddb_hdr_compare,ddb_hdr_free,ddb_hdr_open_read
!!      ddb_hdr_open_write,ddb_malloc,ddb_write_block,read_block,wrtout
!!
!! SOURCE

subroutine mblktyp5 (chkopt,ddbun,dscrpt,filnam,mddb,msym,nddb,vrsddb)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: ddbun,mddb,nddb,vrsddb
 integer,intent(out) :: msym
 character(len=fnlen),intent(in) :: dscrpt,filnam(mddb+1)

!Local variables -------------------------
!scalars
!Define input and output unit numbers:
 integer,parameter :: ddbuntmp=3
 integer :: chkopt,choice,dimekb,iblok,iblok1,iblok2
 integer :: iddb,ii,lmnmax,matom
 integer :: mband,mblktyp,mblok,mkpt,mpert,msize,mtypat
 integer :: nblok,nblokt
 integer :: temp,tmerge,usepaw
 integer,allocatable :: mgblok(:)
 real(dp),parameter :: qtol=2.0d-8
 real(dp) :: diff
 type(ddb_type) :: ddb
 type(ddb_hdr_type) :: ddb_hdr, ddb_hdr8
!arrays
 real(dp),allocatable :: blkval2(:,:,:,:),kpnt(:,:,:)
 character(len=500) :: message

! *********************************************************************

 ! Make sure there is more than one ddb to be read
 if(nddb==1)then
   write(message, '(a,a,a,a,a)' )&
&   'The initialisation mode of MRGDDB, that uses nddb=1,',&
&   'has been disabled in version 2.1 of ABINIT.',&
&   'Action: you should use DDBs that include the symmetry',&
&   'information (and that can be used and merged without',&
&   'initialisation), or you should use ABINITv2.0.'
   MSG_ERROR(message)
 end if

!Evaluate the maximal dimensions of arrays
 dimekb=0 ; matom=0 ; mband=0  ; mblok=0 ; mkpt=0
 msize=0  ; mtypat=0 ; lmnmax=0 ; usepaw=0 ; mblktyp = 1
 msym=192

 do iddb=1,nddb
   call ddb_hdr_open_read(ddb_hdr, filnam(iddb+1), ddbun, vrsddb, dimonly=1)

   mblok=mblok+ddb_hdr%nblok
   mblktyp=max(mblktyp,ddb_hdr%mblktyp)
   matom=max(matom,ddb_hdr%matom)
   mkpt=max(mkpt,ddb_hdr%mkpt)
   mtypat=max(mtypat,ddb_hdr%mtypat)
   msym=max(msym,ddb_hdr%msym)
   mband=max(mband,ddb_hdr%mband)
   dimekb=max(dimekb,ddb_hdr%psps%dimekb)
   lmnmax=max(lmnmax,ddb_hdr%psps%lmnmax)
   usepaw=max(usepaw,ddb_hdr%usepaw)


   call ddb_hdr_free(ddb_hdr)

 end do

 mpert=matom+MPERT_MAX
 msize=3*mpert*3*mpert
 if(mblktyp==3)msize=msize*3*mpert

!write(std_out,*),'msize',msize,'mpert',mpert,'mblktyp',mblktyp
 call ddb%malloc(msize, mblok, matom, mtypat)

!Allocate arrays
 ABI_ALLOCATE(mgblok,(mblok))


!**********************************************************************

!Read the first database

 write(std_out,*)' read the input derivative database information'
 call ddb_hdr_open_read(ddb_hdr, filnam(2), ddbun, vrsddb, &
& matom=matom,mtypat=mtypat,mband=mband,mkpt=mkpt,&
& msym=msym,dimekb=dimekb,lmnmax=lmnmax,usepaw=usepaw)

 ABI_ALLOCATE(blkval2,(2,msize,ddb_hdr%nband(1),mkpt))
 ABI_ALLOCATE(kpnt,(3,mkpt,mblok))

 nblok = ddb_hdr%nblok

 if(nblok>=1)then
!  Read the blocks from the input database.
   write(message, '(a,i5,a)' ) ' read ',nblok,' blocks from the input DDB '
   call wrtout(std_out,message,'COLL')
   choice=1
   do iblok=1,nblok
     call ddb%read_block(iblok,ddb_hdr%nband(1),mpert,&
&     msize,ddb_hdr%nkpt,ddbun,blkval2(1,1,1,1),kpnt(1,1,iblok))
!    Setup merged indicator
     mgblok(iblok)=1
   end do
 else
   call wrtout(std_out,' No bloks in the first ddb ','COLL')
 end if
!Close the first ddb
 close(ddbun)

!*********************************************

!In case of merging of DDBs, iterate the reading
 do iddb=2,nddb

!  Open the corresponding input DDB,
!  and read the database file information
   write(message, '(a,a,i6)' )ch10,&
&   ' read the input derivative database number',iddb
   call wrtout(std_out,message,'COLL')

   call ddb_hdr_open_read(ddb_hdr8, filnam(iddb+1), ddbun, vrsddb, &
&   matom=matom,mtypat=mtypat,mband=mband,mkpt=mkpt,&
&   msym=msym,dimekb=dimekb,lmnmax=lmnmax,usepaw=usepaw)

   if (chkopt==1)then
!    Compare the current DDB and input DDB information.
!    In case of an inconsistency, halt the execution.
     write(message, '(a)' )' compare the current and input DDB information'
     call wrtout(std_out,message,'COLL')

     call ddb_hdr_compare(ddb_hdr, ddb_hdr8)

   else if(chkopt==0)then
!    No comparison between the current DDB and input DDB information.
     write(message, '(a)' )' no comparison between the current and input DDB information'
     call wrtout(std_out,message,'COLL')
     write(message, '(a,a,a)' )&
&     'No comparison/check is performed for the current and input DDB information ',&
&     'because argument --nostrict was passed to the command line. ',&
&     'Use at your own risk !'
     MSG_COMMENT(message)
   end if

   call wrtout(std_out,' Will try to merge this input DDB with the current one.','COLL')

!  First estimate of the total number of bloks, and error
!  message if too large
   write(message, '(a,i5)' ) ' Current number of bloks =',nblok
   call wrtout(std_out,message,'COLL')
   write(message, '(a,i5,a)' )' Will read ',ddb_hdr8%nblok,' blocks from the input DDB '
   call wrtout(std_out,message,'COLL')
   nblokt=nblok+ddb_hdr8%nblok
   if(nblokt>mblok)then
     write(message, '(a,i0,a,a,a,i0,a)' )&
&     'The expected number of blocks',nblokt,' is larger than',ch10,&
&     'the maximum number of blocks',mblok,'.'
     MSG_ERROR(message)
   end if

!  Read the bloks from the temporary database, and close it.
!  Also setup the merging indicator
   choice=1
   do iblok=nblok+1,nblokt
     call ddb%read_block(iblok,ddb_hdr8%nband(1),mpert,&
&     msize,ddb_hdr8%nkpt,ddbun,blkval2(1,1,1,1),kpnt(1,1,iblok))
     mgblok(iblok)=1
   end do
   close(ddbun)

   nblok=nblokt
   write(message, '(a,i5)' ) ' Now, current number of bloks =',nblok
   call wrtout(std_out,message,'COLL')

   ! In certain cases, the different DDB will have different information
   ! on the pseudos (depending on fullinit)
   ! Here, we copy the information of the last DDB file,
   ! only to make the tests pass...
   ddb_hdr%psps%indlmn(:,:,:) = ddb_hdr8%psps%indlmn(:,:,:)
   ddb_hdr%psps%pspso(:) = ddb_hdr8%psps%pspso(:)
   ddb_hdr%psps%ekb(:,:) = ddb_hdr8%psps%ekb(:,:)

   call ddb_hdr_free(ddb_hdr8)

 end do

 call wrtout(std_out,' All DDBs have been read ','COLL')

!*********************************************************

!Check the equality of blocks, and eventually merge them

 if(nblok>=1)then
   call wrtout(std_out,' check the equality of blocks, and eventually merge ','COLL')
   do iblok2=2,nblok
     do iblok1=1,iblok2-1
!      Check the block type identity
       if(ddb%typ(iblok1)==ddb%typ(iblok2))then
!        Check the wavevector identities
         diff=abs(ddb%qpt(1,iblok1)-ddb%qpt(1,iblok2))
         diff=diff+abs(ddb%qpt(2,iblok1)-ddb%qpt(2,iblok2))
         diff=diff+abs(ddb%qpt(3,iblok1)-ddb%qpt(3,iblok2))
         if(abs(diff)<qtol)mgblok(iblok2)=0
       end if
     end do
   end do

!  Count the final number of bloks
   tmerge=0
   do ii=1,nblok
     if(mgblok(ii)==1)tmerge=tmerge+1
   end do
   temp = nblok-tmerge
   nblok=tmerge

!  Summarize the merging phase
   write(message, '(i6,a,i6,a)' )&
&   temp,' blocks are merged; the new DDB will have ',nblok,' blocks.'
   call wrtout(std_out,message,'COLL')
 end if

!**********************************************************************

 write(message, '(a,a)' )' open the output database, write the',' preliminary information '
 call wrtout(std_out,message,'COLL')

 ddb_hdr%dscrpt = trim(dscrpt)
 ddb_hdr%nblok = nblok !nblokt
 ddb_hdr%mblktyp = mblktyp

 call ddb_hdr_open_write(ddb_hdr, filnam(1), ddbun, fullinit=1)

 if(nddb>1)then

!  Write the whole database
   call wrtout(std_out,' write the DDB ','COLL')
   ii = 1 !unit indicator of what will be merged
!  Create a temporary file to decrease memory need.
   do iddb=1,nddb
     call ddb_hdr_open_read(ddb_hdr8, filnam(iddb+1), ddbuntmp, vrsddb)

     do iblok=1,ddb_hdr8%nblok
       if(mgblok(ii)==1) then
         call ddb%read_block(ii,ddb_hdr8%nband(1),mpert,&
&         msize,ddb_hdr8%nkpt,ddbuntmp,blkval2(:,:,:,:),kpnt(:,:,ii))
         choice=2
         call ddb%write_block(ii,choice,ddb_hdr%nband(1),mpert,&
&         msize,ddb_hdr8%nkpt,ddbun,blkval2(:,:,:,:),kpnt(:,:,ii))
       else
         write(message, '(a,i4,a,i4,a)' )&
&         ' Bloc number',iblok,' of DDB ',iddb,&
&         ' was merged, so do not write it'
         call wrtout(std_out,message,'COLL')
       end if
       ii = ii+1
     end do
     close(ddbuntmp)
     call ddb_hdr_free(ddb_hdr8)
   end do !iddb=1,nddb

!  Also write summary of bloks at the end
   write(ddbun, '(/,a)' )' List of bloks and their characteristics '
   choice=3
   do iblok=1,nblokt
     if(mgblok(iblok)==1)then
       call ddb%write_block(iblok,choice,ddb_hdr%nband(1),mpert,&
&       msize,ddb_hdr%nkpt,ddbun)
     end if
   end do

 end if

 close (ddbun)

!*********************************************************************

 call ddb_hdr_free(ddb_hdr)

!Deallocate arrays
 ABI_DEALLOCATE(mgblok)
 ABI_DEALLOCATE(blkval2)
 ABI_DEALLOCATE(kpnt)

 call ddb_free(ddb)

end subroutine mblktyp5
!!***

#ifdef MR_DEV
!!****f* m_ddb/dfpt_lw_doutput
!! NAME
!! dfpt_lw_doutput
!!
!! FUNCTION
!! Write the matrix of third-order derivatives from the long wave calculation 
!! to the to the DDB
!!
!! INPUTS
!!  blkflg(3,mpert,3,mpert,3,mpert)= ( 1 if the element of the 3dte
!!   has been calculated ; 0 otherwise )
!!  d3(2,3,mpert,6,mpert,3,mpert)= matrix of the 3DTE
!!  mpert =maximum number of ipert
!!  natom=Number of atoms
!!  ntypat=Number of type of atoms
!!  unddb = unit number for DDB output
!!
!! NOTES
!! - d3 holds the third-order derivatives before computing
!!   the permutations of the perturbations.
!! - the dimension 6 of the 4th argument in d3 is used to define
!!   both up and down extradiagonal strain perturbations.                          
!!   Necessary because their q-gradient is not symmetric. 
!!
!! PARENTS
!!      longwave
!!
!! CHILDREN
!!      ddb_free,ddb_malloc,ddb_write_blok,wrtout
!!
!! SOURCE

subroutine dfpt_lw_doutput(blkflg,d3,mpert,natom,ntypat,unddb)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: mpert,unddb,natom,ntypat
!arrays
 integer,intent(in) :: blkflg(3,mpert,3,mpert,3,mpert)
 real(dp),intent(in) :: d3(2,3,mpert,3,mpert,3,mpert)

!Local variables -------------------------
!scalars
 integer :: idir1,idir2,idir3,ii,index,ipert1,ipert2,ipert3,msize,nelmts
 character(len=500) :: message
 type(ddb_type) :: ddb

!*************************************************************************

 msize = 54*mpert*mpert*mpert
 call ddb_malloc(ddb,msize,1,natom,ntypat)

 ddb%nrm = one
 ddb%qpt = zero  

 index=0
 do ipert3=1,mpert
   do idir3=1,3
     do ipert2=1,mpert
       do idir2=1,3
         do ipert1=1,mpert
           do idir1=1,3
             index=index + 1
             ddb%flg(index,1) = blkflg(idir1,ipert1,idir2,ipert2,idir3,ipert3)
             ddb%val(:,index,1)= d3(:,idir1,ipert1,idir2,ipert2,idir3,ipert3)

           end do
         end do
       end do
     end do
   end do
 end do

!########  Write blok of third-order derivatives to DDB

!Count the number of elements
 nelmts=0
 do ii=1,msize
   if(ddb%flg(ii,1)==1)nelmts=nelmts+1
 end do

!Write the block type and number of elements
 write(unddb,*)' '
 write(unddb, '(a,i8)' )' 3rd derivatives              - # elements :',nelmts

!Write the phonon wavevectors
 write(unddb, '(a,3es16.8,f6.1)' )' qpt',(ddb%qpt(ii,1),ii=1,3),ddb%nrm(1,1)
 write(unddb, '(a,3es16.8,f6.1)' )'    ',(ddb%qpt(ii,1),ii=4,6),ddb%nrm(2,1)
 write(unddb, '(a,3es16.8,f6.1)' )'    ',(ddb%qpt(ii,1),ii=7,9),ddb%nrm(3,1)

!Write the matrix elements 
 index=0
 do ipert3=1,mpert
   do idir3=1,3
     do ipert2=1,mpert
       do idir2=1,3
         do ipert1=1,mpert
           do idir1=1,3
             index=index+1
             if (ddb%flg(index,1)==1)then
               write(unddb, '(6i4,2d22.14)' )&
&              idir1,ipert1,idir2,ipert2,idir3,ipert3,ddb%val(1,index,1),ddb%val(2,index,1)
             end if
           end do
         end do
       end do
     end do
   end do
 end do

end subroutine dfpt_lw_doutput
!!***
#endif

end module m_ddb
!!***

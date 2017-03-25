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
!! Copyright (C) 2011-2017 ABINIT group (MJV, XG, MT, MM, MVeithen, MG, PB, JCC)
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
 use m_profiling_abi
 use m_errors
 use m_xmpi

 use m_fstrings,       only : sjoin, itoa, ktoa
 use m_numeric_tools,  only : mkherm
 use m_io_tools,       only : open_file, get_unit
 use m_copy,           only : alloc_copy
 use m_geometry,       only : phdispl_cart2red
 use m_crystal,        only : crystal_t, crystal_init
 use m_pawtab,         only : pawtab_type,pawtab_nullify,pawtab_free
 use m_dynmat,         only : cart29, d2sym3, cart39, d3sym, chneu9, asria_calc, asria_corr, asrprs, dfpt_phfrq

 implicit none

 private

 public :: ddb_getdims      ! Open a DDB file and read basic dimensions and variables.
 public :: gtblk9           ! Finds the block containing the derivatives of the total energy.
 public :: read_blok8       ! This routine reads blocks of data in the DDBs.
 public :: psddb8           ! I/O of pseudopotentials. Read also the number of data blocks.
 public :: ioddb8_in        ! Open Derivative DataBase, and read preliminary information.
 public :: rdddb9           ! This routine reads the derivative database entirely,
 public :: nlopt            ! Output of all quantities related to third-order derivatives of the energy.
 public :: chkin9
 public :: carttransf       ! Transform a second-derivative matrix (EIG2D) from reduced
                            ! coordinates to cartesian coordinates.

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

 end type ddb_type

 public :: ddb_from_file            ! Construct the object from the DDB file.
 public :: ddb_free                 ! Free dynamic memory.
 public :: ddb_malloc               ! Allocate dynamic memory
 public :: ddb_bcast                ! Broadcast the object.
 public :: ddb_copy                 ! Copy the object.
 public :: ddb_get_etotal           ! Read the GS total energy.
 public :: ddb_get_dielt_zeff       ! Reads the Dielectric Tensor and the Effective Charges
 public :: ddb_get_dchidet          ! Reads the non-linear optical susceptibility tensor and the
                                    ! first-order change in the linear dielectric susceptibility
 public :: ddb_diagoq               ! Compute the phonon frequencies at the specified q-point by performing
                                    ! a direct diagonalizatin of the dynamical matrix.
 public :: ddb_get_asrq0            ! Return object used to enforce the acoustic sum rule
                                    ! from the Dynamical matrix at Gamma. Used in ddb_diagoq.

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

 end type asrq0_t

 public :: asrq0_apply      ! Impose the acoustic sum rule based on the q=0 block found in the DDB file.
 public :: asrq0_free       ! Free memory
!!***

 ! TODO: We should use this constants instead of magic numbers!
 ! BTW: Using a different value for NOSTAT and STAT is a non-sense!
 ! They both are 2-th order derivatives of the total energy!

 ! Flags used to indentify the block type.
 !integer,private,parameter :: BLKTYPE_ETOT = 0         ! Total energy
 !integer,private,parameter :: BLKTYPE_2DE_NOSTAT = 1   ! Second order derivative of the energy (non-stationary expression)
 !integer,private,parameter :: BLKTYPE_2DE_STAT = 2     ! Second order derivative of the energy (stationary expression)
 !integer,private,parameter :: BLKTYPE_3DE = 3          ! Third order derivative of the energy
 !integer,private,parameter :: BLKTYPE_1DE = 4          ! First order derivative of the energy
 !integer,private,parameter :: BLKTYPE_2DEIG = 5        ! Second order derivative of the eigenvalues

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
!!      anaddb,dfpt_looppert,dfptnl_doutput,eph,gstate
!!      m_effective_potential_file,mblktyp1,mblktyp5,thmeig
!!
!! CHILDREN
!!
!! SOURCE

subroutine ddb_free(ddb)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ddb_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(ddb_type),intent(inout) :: ddb

! ************************************************************************

 !integer
 if (allocated(ddb%flg))  then
   ABI_FREE(ddb%flg)
 end if
 if (allocated(ddb%typ))  then
   ABI_FREE(ddb%typ)
 end if

 ! real
 if (allocated(ddb%amu))  then
   ABI_FREE(ddb%amu)
 end if
 if (allocated(ddb%nrm))  then
   ABI_FREE(ddb%nrm)
 end if
 if (allocated(ddb%qpt))  then
   ABI_FREE(ddb%qpt)
 end if
 if (allocated(ddb%val))  then
   ABI_FREE(ddb%val)
 end if

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
!!
!! SOURCE

subroutine ddb_copy(iddb, oddb)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ddb_copy'
!End of the abilint section

 implicit none

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
 call alloc_copy(iddb%flg,oddb%flg)
 call alloc_copy(iddb%typ,oddb%typ)
 call alloc_copy(iddb%amu,oddb%amu)
 call alloc_copy(iddb%nrm,oddb%nrm)
 call alloc_copy(iddb%qpt,oddb%qpt)
 call alloc_copy(iddb%val,oddb%val)

end subroutine ddb_copy
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_malloc
!! NAME
!! ddb_malloc
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      dfptnl_doutput,gstate,m_ddb,mblktyp1,mblktyp5,thmeig
!!
!! CHILDREN
!!
!! SOURCE

subroutine ddb_malloc(ddb,msize,nblok,natom,ntypat)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ddb_malloc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!array
 integer,intent(in) :: msize,nblok,natom,ntypat
 type(ddb_type),intent(inout) :: ddb

! ************************************************************************

 ! FIXME
 ! This is done in rdddb9 but only by the master node!
 ! Should rationalize the different steps
 ddb%msize = msize
 ddb%nblok = nblok
 ddb%natom = natom
 ddb%mpert = natom+6
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
!!
!! SOURCE

subroutine ddb_bcast(Ddb, master, comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ddb_bcast'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!array
 type(ddb_type),intent(inout) :: ddb
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
   call ddb_malloc(ddb,ddb%msize,ddb%nblok,ddb%natom,ddb%ntypat)
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

!!****f* m_ddb/inprep8
!! NAME
!! inprep8
!!
!! FUNCTION
!! Open Derivative DataBase, then reads the variables that
!! must be known in order to dimension the arrays before complete reading
!! Note: only one processor read or write the DDB.
!!
!! INPUTS
!! character(len=*) filnam: name of input or output file
!! unddb=unit number for input or output
!! vrsddb=6 digit integer giving date, in form yymmdd for month=mm(1-12),
!!  day=dd(1-31), and year=yy(90-99 for 1990 to 1999,00-89 for 2000 to 2089), of current DDB version.
!!
!! OUTPUT
!! dimekb=dimension of ekb (only used for norm-conserving psps)
!! lmnmax=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!       =if useylm=0, max number of (l,n)   comp. over all type of psps
!! mband=maximum number of bands
!! mblktyp=largest block type
!! msym=maximum number of symmetries
!! natom=number of atoms
!! nblok=number of bloks in the DDB
!! nkpt=number of k points
!! ntypat=number of atom types
!! usepaw= 0 for non paw calculation; =1 for paw calculation
!!
!! PARENTS
!!      m_ddb
!!
!! CHILDREN
!!
!! SOURCE


subroutine inprep8 (dimekb,filnam,lmnmax,mband,mblktyp,msym,natom,nblok,nkpt,&
& ntypat,unddb,usepaw,vrsddb)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'inprep8'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: unddb,vrsddb
 integer, intent(out) :: msym
 integer,intent(out) :: dimekb,lmnmax,mband,mblktyp,natom,nblok,nkpt,ntypat,usepaw
 character(len=*),intent(in) :: filnam

!Local variables -------------------------
!scalars
!Set routine version number here:
 integer,parameter :: vrsio8=100401,vrsio8_old=010929,vrsio8_old_old=990527
 integer :: bantot,basis_size0,blktyp,ddbvrs,iband,iblok,iekb,ii,ikpt,iline,im,ios,iproj
 integer :: itypat,itypat0,jekb,lmn_size0,mproj,mpsang,nekb,nelmts,nsppol
 integer :: occopt,pspso0,nsym
 logical :: ddbvrs_is_current_or_old,testn,testv
 character(len=12) :: string
 character(len=32) :: blkname
 character(len=500) :: message
 character(len=6) :: name_old
 character(len=80) :: rdstring
!arrays
 integer,allocatable :: nband(:)
 character(len=12) :: name(9)

! *********************************************************************

!Check inprep8 version number (vrsio8) against mkddb version number (vrsddb)
 if (vrsio8/=vrsddb) then
   write(message, '(a,i0,2a,i0)' )&
&   'The input/output DDB version number= ',vrsio8,ch10,&
&   'is not equal to the DDB version number= ',vrsddb
   MSG_BUG(message)
 end if

!Open the input derivative database.
 call wrtout(std_out, sjoin(" Opening DDB file:", filnam), 'COLL')
 if (open_file(filnam,message,unit=unddb,form="formatted",status="old",action="read") /= 0) then
   MSG_ERROR(message)
 end if

!Check the compatibility of the input DDB with the DDB code
 read (unddb,*)
 read (unddb,*)
 read (unddb, '(20x,i10)' )ddbvrs

 if (all(ddbvrs/= [vrsio8, vrsio8_old, vrsio8_old_old]) )then
   write(message, '(a,i10,2a,3(a,i10))' )&
&   'The input DDB version number=',ddbvrs,' does not agree',ch10,&
&   'with the allowed code DDB version numbers,',vrsio8,', ',vrsio8_old,' and ',vrsio8_old_old
   MSG_ERROR(message)
 end if

!Read the 4 n-integers, also testing the names of data,
!and checking that their value is acceptable.
!This is important to insure that any array has a sufficient dimension.
 read (unddb,*)
 read (unddb,*)
 read (unddb,*)
 testn=.true.
 testv=.true.
 ddbvrs_is_current_or_old=(ddbvrs==vrsio8.or.ddbvrs==vrsio8_old)

!1. usepaw
 if(ddbvrs==vrsio8)then
   read (unddb, '(1x,a9,i10)' )name(1),usepaw
 else
   usepaw=0;name(1)='   usepaw'
 end if
 if(name(1)/='   usepaw')testn=.false.
!2. natom
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,i10)' )name(2),natom
 else
   read (unddb, '(1x,a6,i10)' )name_old,natom ; name(2)='   '//name_old
 end if
 if(name(2)/='    natom')testn=.false.
 if(natom<=0)testv=.false.
!3. nkpt
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,i10)' )name(3),nkpt
 else
   read (unddb, '(1x,a6,i10)' )name_old,nkpt ; name(3)='   '//name_old
 end if
 if(name(3)/='     nkpt')testn=.false.
 if(nkpt <=0)testv=.false.
!4. nsppol
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,i10)' )name(4),nsppol
 else
   read (unddb, '(1x,a6,i10)' )name_old,nsppol ; name(4)='   '//name_old
 end if
 if(name(4)/='   nsppol')testn=.false.
 if(nsppol<=0.or.nsppol>2)testv=.false.
!5. nsym
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,i10)' )name(5),nsym
 else
   read (unddb, '(1x,a6,i10)' )name_old,nsym ; name(5)='   '//name_old
 end if
 if(name(5)/='     nsym')testn=.false.
!MG FIXME Why this and why do we need msym?
 msym = 192
 if (nsym > msym) msym=nsym
!if(nsym <=0.or.nsym >msym )testv=.false.
!6. ntypat
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,i10)' )name(6),ntypat
 else
   read (unddb, '(1x,a6,i10)' )name_old,ntypat ; name(6)='   '//name_old
 end if
 if(name(6)/='   ntypat' .and. name(6)/='    ntype')testn=.false.
 if(ntypat<=0)testv=.false.
!7. occopt
!Before reading nband, the last parameters that define
!the dimension of some array, need to know what is their
!representation, given by occopt
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,i10)' )name(7),occopt
 else
   read (unddb, '(1x,a6,i10)' )name_old,occopt ; name(7)='   '//name_old
 end if
 if(name(7)/='   occopt')testn=.false.
 if(occopt<0.or.occopt>8)testv=.false.

!Message if the names or values are not right
 if (.not.testn.or..not.testv) then
   write(message, '(a,a)' )' inprep8 : An error has been found in the',' positive n-integers contained in the DDB : '
   call wrtout(std_out,message,'COLL')
   write(message, '(a)' )   '     Expected                      Found     '
   call wrtout(std_out,message,'COLL')
   write(message, '(a,i9,a,a,a,i10)' )'    natom , larger than',0,'    ',trim(name(2)),' =',natom
   call wrtout(std_out,message,'COLL')
   write(message, '(a,i9,a,a,a,i10)' )'    nkpt  , larger than',0,'    ',trim(name(3)),' =',nkpt
   call wrtout(std_out,message,'COLL')
   write(message, '(a,i1,a,a,a,i10)' )'    nsppol, either    1 or     ',2,'    ',trim(name(4)),' =',nsppol
   call wrtout(std_out,message,'COLL')
!  write(message, '(a,i10,a,a,a,i10)' )&   '    nsym  , lower than',msym,'    ',trim(name(5)),' =',nsym
   call wrtout(std_out,message,'COLL')
   write(message, '(a,i9,a,a,a,i10)' )'    ntypat , larger than',0,'   ',trim(name(6)),' =',ntypat
   call wrtout(std_out,message,'COLL')
   write(message, '(a,a,a,i10)' )'    occopt,     equal to 0,1 or 2   ',trim(name(7)),' =',occopt
   call wrtout(std_out,message,'COLL')

   MSG_ERROR('See the error message above.')
 end if

!One more set of parameters define the dimensions of the
!array : nband. Morever, it depends on occopt and !nkpt, and has to be
!tested after the test on nkpt is performed.

!8. nband
 ABI_MALLOC(nband,(nkpt))
 if(occopt==2)then
   im=12
   do iline=1,(nkpt+11)/12
     if(iline==(nkpt+11)/12)im=nkpt-12*(iline-1)
     if(ddbvrs_is_current_or_old)then
       read (unddb, '(1x,a9,5x,12i5)' )name(1),(nband((iline-1)*12+ii),ii=1,im)
     else
       read (unddb, '(1x,a6,5x,12i5)' )name_old,&
&       (nband((iline-1)*12+ii),ii=1,im) ; name(1)='   '//name_old
     end if
     if (iline==1) then
       call ddb_chkname(name(1),'    nband')
     else
       call ddb_chkname(name(1),'         ')
     end if
   end do
 else
   if(ddbvrs_is_current_or_old)then
     read (unddb, '(1x,a9,i10)' )name(1),nband(1)
   else
     read (unddb, '(1x,a6,i10)' )name_old,nband(1) ; name(1)='   '//name_old
   end if
   call ddb_chkname(name(1),'    nband')
   if(nkpt>1)then
     do ikpt=2,nkpt
       nband(ikpt)=nband(1)
     end do
   end if
 end if

!Check all nband values, and sum them
 bantot=0
 do ikpt=1,nkpt
   if(nband(ikpt)<0)then
     write(message, '(a,i0,a,i0,3a)' )&
&     'For ikpt = ',ikpt,'  nband = ',nband(ikpt),' is negative.',ch10,&
&     'Action: correct your DDB.'
     MSG_ERROR(message)
   end if
   bantot=bantot+nband(ikpt)
 end do

 mband=maxval(nband(:))

!Skip the rest of variables
!9. acell
 read (unddb,*)
!10. amu
 do iline=1,(ntypat+2)/3
   read (unddb,*)
 end do
!11. dilatmx
 if(ddbvrs_is_current_or_old) read (unddb,*)
!12. ecut
 read (unddb,*)
!12b. pawecutdg (PAW only)
 if(ddbvrs==vrsio8.and.usepaw==1) read (unddb,*)
!13. ecutsm
 if(ddbvrs_is_current_or_old) read (unddb,*)
!14. intxc
 if(ddbvrs_is_current_or_old) read (unddb,*)
!15. iscf
 read (unddb,*)
!16. ixc
 read (unddb,*)
!17. kpt
 do iline=1,nkpt
   read (unddb,*)
 end do
!18. kptnrm
 read (unddb,*)
!19. ngfft
 read (unddb,*)
!20. nspden
 if(ddbvrs_is_current_or_old) read (unddb,*)
!21. nspinor
 if(ddbvrs_is_current_or_old) read (unddb,*)
!22. occ
 if(occopt==2)then
   do iline=1,(bantot+2)/3
     read (unddb,*)
   end do
 else
   write(message,*)' inprep8 : nband(1)=',nband(1)
   call wrtout(std_out,message,'COLL')
   do iline=1,(nband(1)+2)/3
     read (unddb,'(a80)')rdstring
     write(message,*)trim(rdstring)
     call wrtout(std_out,message,'COLL')
   end do
 end if
!23. rprim
 do iline=1,3
   read (unddb,*)
 end do
!24. dfpt_sciss
 read (unddb,*)
!25. spinat
 if(ddbvrs_is_current_or_old)then
   do iline=1,natom
     read (unddb,*)
   end do
 end if
!26. symafm
 if(ddbvrs_is_current_or_old)then
   do iline=1,(nsym+11)/12
     read (unddb,*)
   end do
 end if
!27. symrel
 do iline=1,nsym
   read (unddb,*)
 end do
!28old. xred
 if(.not.ddbvrs_is_current_or_old)then
   do iline=1,natom
     read (unddb,*)
   end do
 end if
!28. tnons
 do iline=1,nsym
   read (unddb,*)
 end do
!29. tolwfr
 if(ddbvrs_is_current_or_old) read (unddb,*)
!30. tphysel
 if(ddbvrs_is_current_or_old) read (unddb,*)
!31. tsmear
 if(ddbvrs_is_current_or_old) read (unddb,*)
!32. type
 do iline=1,(natom+11)/12
   read (unddb,*)
 end do
!33old. tolwfr
 if(.not.ddbvrs_is_current_or_old) read (unddb,*)
!33. wtk
 do iline=1,(nkpt+2)/3
   read (unddb,*)
 end do
!34. xred
 if(ddbvrs_is_current_or_old)then
   do iline=1,natom
     read (unddb,*)
   end do
 end if
!35. znucl
 if(ddbvrs_is_current_or_old)then
   do iline=1,(ntypat+2)/3
     read (unddb,*)
   end do
 end if
!36. zion
 do iline=1,(ntypat+2)/3
   read (unddb,*)
 end do

 read (unddb,*)

!Now, take care of the pseudopotentials
 read(unddb, '(a12)' )string

 if(string=='  Descriptio')then

   read (unddb,*)
   if (ddbvrs==vrsio8_old.or.ddbvrs==vrsio8_old_old) then
     read (unddb, '(10x,i3,14x,i3,11x,i3)', iostat=ios )dimekb,lmnmax,usepaw
     if(ios/=0)then
       backspace(unddb)
       read (unddb, '(10x,i3,14x,i3)')dimekb,lmnmax
       usepaw=0
     end if
   else if (ddbvrs==vrsio8) then
     read (unddb, '(10x,i3)') usepaw
     if (usepaw==0) then
       read (unddb, '(10x,i3,14x,i3)' ) dimekb,lmnmax
     else
       dimekb=0;lmnmax=0
     end if
   end if
   if (usepaw==0) then
     do itypat=1,ntypat
       read(unddb, '(13x,i4,9x,i3,8x,i4)' )itypat0,pspso0,nekb
       read(unddb,*)
       do iekb=1,nekb
         do jekb=1,nekb,4
           read(unddb,*)
         end do
       end do
     end do
   else
     do itypat=1,ntypat
       read(unddb, '(12x,i4,12x,i3,12x,i5)' )itypat0,basis_size0,lmn_size0
       lmnmax=max(lmnmax,lmn_size0)
       read(unddb,*)
       read(unddb,*)
       read(unddb,'(24x,i3)') nekb
       read(unddb,*)
       do iekb=1,nekb,4
         read(unddb,*)
       end do
     end do
   end if

 else if(string==' Description')then
   if (usepaw==1) then
     MSG_BUG('old DDB pspformat not compatible with PAW 1')
   end if

   read (unddb, '(10x,i3,10x,i3)' )mproj,mpsang
   dimekb=mproj*mpsang
   usepaw=0
   do itypat=1,ntypat
     read (unddb,*)
!    For f-electrons, one more line has been written
     do iproj=1,mproj*max(1,(mpsang+2)/3)
       read (unddb,*)
     end do
   end do

 else if(string==' No informat')then

   dimekb=0
   lmnmax=0
   usepaw=0

 else
   write(message, '(a,a,a,a)' )&
&   'Error when reading the psp information',ch10,&
&   'String=',trim(string)
   MSG_BUG(message)
 end if

!Now, the number of blocks
 read(unddb,*)
 read(unddb,*)
 read(unddb, '(24x,i4)' )nblok

!Now, the type of each blok, in turn
 mblktyp=1
 if(nblok>=1)then
   do iblok=1,nblok

     read(unddb,*)
     read(unddb, '(a32,12x,i8)' )blkname,nelmts
     if(blkname==' 2nd derivatives (non-stat.)  - ' .or.  blkname==' 2rd derivatives (non-stat.)  - ')then
       blktyp=1
     else if(blkname==' 2nd derivatives (stationary) - ' .or. blkname==' 2rd derivatives (stationary) - ')then
       blktyp=2
     else if(blkname==' 3rd derivatives              - ')then
       blktyp=3
     else if(blkname==' Total energy                 - ')then
       blktyp=0
     else if(blkname==' 1st derivatives              - ')then
       blktyp=4
     else if(blkname==' 2nd eigenvalue derivatives   - ' .or. blkname==' 2rd eigenvalue derivatives   - ')then
       blktyp=5
     else
       write(message, '(a,a,a,a,a,a,a,a,a)' )&
&       'The following string appears in the DDB in place of',' the block type description :',ch10,blkname,ch10,&
&       'Action: check your DDB.',ch10,&
&       'Note: If you did use an abinit version prior to 6.12 to generate your DDB',&
&       'pay attention to the change:: 2rd derivatives ==> 2nd derivatives'
       MSG_ERROR(message)
     end if

     if(blktyp==1.or.blktyp==2)then
!      Read the phonon wavevector
       read(unddb,*)
     else if(blktyp==3)then
!      Read the perturbation wavevectors
       read(unddb,*)
       read(unddb,*)
       read(unddb,*)
       mblktyp=3
     else if(blktyp==5)then
       read(unddb,*)
       mblktyp=5
     end if

!    Read every element
     if(blktyp==5)then
       do ikpt=1,nkpt
         read(unddb,*)
         do iband=1,nband(ikpt)
           read(unddb,*)
           do ii=1,nelmts
             read(unddb,*)
           end do
         end do
       end do
     else
       do ii=1,nelmts
         read(unddb,*)
       end do
     end if

   end do
 end if

 ABI_FREE(nband)

!Close the DDB
 close(unddb)

end subroutine inprep8
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ddb_getdims
!! NAME
!! ddb_getdims
!!
!! FUNCTION
!! Open Derivative DataBase, then reads the variables that
!! must be known in order to dimension the arrays before complete reading
!!
!! INPUTS
!! character(len=*) filnam: name of input or output file
!! unddb=unit number for input or output
!! vrsddb=6 digit integer giving date, in form yymmdd for month=mm(1-12),
!!  day=dd(1-31), and year=yy(90-99 for 1990 to 1999,00-89 for 2000 to 2089),
!!  of current DDB version.
!!
!! OUTPUT
!! dimekb=dimension of ekb (only used for norm-conserving psps)
!! lmnmax=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!       =if useylm=0, max number of (l,n)   comp. over all type of psps
!! mband=maximum number of bands
!! mblktyp=largest block type
!! msym=maximum number of symmetries
!! natom=number of atoms
!! nblok=number of bloks in the DDB
!! nkpt=number of k points
!! ntypat=number of atom types
!! usepaw= 0 for non paw calculation; =1 for paw calculation
!! comm=MPI communicator.
!!
!! PARENTS
!!      anaddb,m_ddb,m_effective_potential_file,mblktyp1,mblktyp5,mrgddb
!!
!! CHILDREN
!!
!! SOURCE

subroutine ddb_getdims(dimekb,filnam,lmnmax,mband,mblktyp,msym,natom,nblok,nkpt,ntypat,unddb,usepaw,vrsddb,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ddb_getdims'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: unddb,vrsddb,comm
 integer,intent(out) :: msym,dimekb,lmnmax,mband,mblktyp,natom,nblok,nkpt,ntypat,usepaw
 character(len=*),intent(in) :: filnam

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: ierr
 !integer :: mpert,msize

! *********************************************************************

 ! Master node reads dims from file and then broadcast.
 if (xmpi_comm_rank(comm) == master) then
   call inprep8(dimekb,filnam,lmnmax,mband,mblktyp,msym,natom,nblok,nkpt,ntypat,unddb,usepaw,vrsddb)
 end if

 if (xmpi_comm_size(comm) > 1) then
   call xmpi_bcast(dimekb, master, comm, ierr)
   call xmpi_bcast(lmnmax, master, comm, ierr)
   call xmpi_bcast(mband, master, comm, ierr)
   call xmpi_bcast(mblktyp, master, comm, ierr)
   call xmpi_bcast(msym, master, comm, ierr)
   call xmpi_bcast(natom, master, comm, ierr)
   call xmpi_bcast(nblok, master, comm, ierr)
   call xmpi_bcast(nkpt, master, comm, ierr)
   call xmpi_bcast(ntypat, master, comm, ierr)
   call xmpi_bcast(usepaw, master, comm, ierr)
 end if

 ! Maximum number of perturbations and size of matrix.
 !mpert=natom+6
 !msize=3*mpert*3*mpert; if (mblktyp==3) msize=msize*3*mpert

end subroutine ddb_getdims
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/gtblk9
!!
!! NAME
!! gtblk9
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
!! rfphon(4) = 1=> response to phonons (for the four possible derivatives. Two should be used for a second derivative of total energy)
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
!!      anaddb,m_ddb,m_effective_potential_file,m_phonons,thmeig
!!
!! CHILDREN
!!
!! SOURCE


subroutine gtblk9(ddb,iblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,rftyp)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gtblk9'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: rftyp
 integer,intent(out) :: iblok
 type(ddb_type),intent(in) :: ddb
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
   !  Identifies if qphon is at gamma
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
           write(message,'(a,i0,a,a,a)' )&
&           'The block ',iblok,' does not match the requirement',ch10,&
&           'because it lacks the total energy'
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
     write(message, '(a,a,a)' )&
&     ' gtblk9 : ',ch10,&
&     '  Unable to find block corresponding to the following specifications :'
     call wrtout(std_out,message,'COLL')
     write(message, '(a,i3)' )' Type (rfmeth) =',rftyp
     call wrtout(std_out,message,'COLL')
     write(message, '(a)' ) ' ider qphon(3)         qphnrm   rfphon rfelfd rfstrs'
     call wrtout(std_out,message,'COLL')
!    write(std_out,*)' nder=',nder
     do ider=1,nder
       write(message, '(i4,4f6.2,3i7)' )&
&       ider,(qphon(ii,ider),ii=1,3),qphnrm(ider),rfphon(ider),rfelfd(ider),rfstrs(ider)
       call wrtout(std_out,message,'COLL')
     end do
   end if
 end if

 if (ok==1 .and. ddb%prtvol > 1) then
   write(message,'(a,i0,a,a)')' gtblk9: found block number ',iblok,' agree with',' specifications '
   call wrtout(std_out,message,'COLL')
 end if

 ABI_FREE(worki)

end subroutine gtblk9
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
!!
!! SOURCE


subroutine gamma9(gamma,qphon,qphnrm,qtol)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gamma9'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(out) :: gamma
 real(dp),intent(in) :: qphnrm,qtol
!arrays
 real(dp),intent(in) :: qphon(3)

! *********************************************************************

 if( (  abs(qphon(1))<qtol .and.  &
& abs(qphon(2))<qtol .and.        &
& abs(qphon(3))<qtol      ) .or.  &
& abs(qphnrm)<qtol )then
   gamma=1
 else
   gamma=0
 end if

end subroutine gamma9
!!***

!----------------------------------------------------------------------

!!****f* m_db_blk/read_blok8
!!
!! NAME
!! read_blok8
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
!!
!! SOURCE

subroutine read_blok8(ddb,iblok,mband,mpert,msize,nkpt,nunit,&
&     blkval2,kpt) !optional


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'read_blok8'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: mband,mpert,msize,nkpt,nunit
 integer, intent(in) :: iblok
 type(ddb_type),intent(inout) :: ddb
!arrays
 real(dp),intent(out),optional :: kpt(3,nkpt)
 real(dp),intent(out),optional :: blkval2(2,msize,mband,nkpt)

!Local variables -------------------------
!scalars
 integer :: band,iband,idir1,idir2,idir3,ii,ikpt,index,ipert1,ipert2,ipert3
 integer :: nelmts
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
   write(message, '(a,a,a,a,a,a)' )&
&   'The following string appears in the DDB in place of',&
&   ' the block type description :',ch10,trim(name),ch10,&
&   'Action: check your DDB.'
   MSG_ERROR(message)
 end if

!Read the 2nd derivative block
 if(ddb%typ(iblok)==1.or.ddb%typ(iblok)==2)then

!  First check if there is enough space to read it
   if(msize<(3*mpert*3*mpert))then
     write(message,'(a,a,a,i10,a,i10,a,a,a)')&
&     'There is not enough space to read a second-derivative block.',ch10,&

&     'Action: increase msize and recompile.'
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
&     'There is not enough space to read a third-derivative block.',ch10,&
&     'The size provided is only ',msize,' although ',3*mpert*3*mpert*3*mpert,' is needed.',ch10,&
&     'Action: increase msize and recompile.'
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
     write(message, '(a,a,a,i10,a,i10,a,a,a)' )&
&     'There is not enough space to read a first-derivative block.',ch10,&
&     'The size provided is only ',msize,' although ',3*mpert,' is needed.',ch10,&
&     'Action: increase msize and recompile.'
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
     write(message, '(a,a,a,i10,a,i10,a,a,a)' )&
&     'There is not enough space to read a second-derivative block.',ch10,&
&     'The size provided is only ',msize,' although ',3*mpert*3*mpert*mband*nkpt,' is needed.',ch10,&
&     'Action: increase msize and recompile.'
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

end subroutine read_blok8
!!***

!----------------------------------------------------------------------

!!****f* m_db_blk/psddb8
!!
!! NAME
!! psddb8
!!
!! FUNCTION
!! Take care of the i/o of pseudopotentials for the
!! Derivative DataBase, and also the number of data blocks.
!!
!! INPUTS
!!  choice=(1 => read), (2=> write)
!!  dimekb=dimension of ekb (contains Kleimann-Bylander energies)
!!         used only for norm-conserving pseudopotentials
!!  lmnmax=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!        =if useylm=0, max number of (l,n)   comp. over all type of psps
!!  nunit=unit number for the Derivative DataBase.
!!  ntypat=number of atom types
!!  pspso(ntypat)=For each type of psp, 1 if no spin-orbit component is taken
!!     into account, 2 if a spin-orbit component is used
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  useylm=governs the way the nonlocal operator is to be applied:
!!         1=using Ylm, 0=using Legendre polynomials
!!  vrsddb=Derivative Database version, for check of compatibility.
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  indlmn(6,i,ntypat)= array giving l,m,n,lm,ln,s for i=ln  (if useylm=0)
!!                                                 or i=lmn  (if useylm=1)
!!  ekb(dimekb,ntypat)= (norm-conserving psps only) (Real) Kleinman-Bylander energies (hartree)
!!                      Presently the only characteristics of the psp
!!  fullinit=0 if the ekb are not available, at input as well as at output
!!  pawtab(ntypat*usepaw)= (PAW only) PAW datasets characteristics
!!                  Presently only pawtab%basis_size,pawtab%lmn_size,pawtab%shape_type
!!                  pawtab%rpaw,pawtab%rshp,pawtab%dij0  are used
!!  nblok=number of blocks
!!
!! NOTES
!! Only executed by one processor
!!
!! PARENTS
!!      dfpt_looppert,eig2tot,gstate,m_ddb,mblktyp1,mblktyp5,nonlinear,respfn
!!      thmeig
!!
!! CHILDREN
!!
!! SOURCE

subroutine psddb8 (choice,dimekb,ekb,fullinit,indlmn,lmnmax,&
&          nblok,ntypat,nunit,pawtab,pspso,usepaw,useylm,vrsddb)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psddb8'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: choice,dimekb,lmnmax,ntypat,nunit,usepaw,useylm
 integer,intent(in) :: vrsddb
 integer,intent(inout) :: fullinit,nblok
!arrays
 integer,intent(in) :: pspso(ntypat)
 integer,intent(inout) :: indlmn(6,lmnmax,ntypat)
 real(dp),intent(inout) :: ekb(dimekb,ntypat)
 type(pawtab_type),intent(inout) :: pawtab(ntypat*usepaw)

!Local variables -------------------------
!Set the version number
!scalars
 integer,parameter :: vrsio8=100401,vrsio8_old=010929,vrsio8_old_old=990527
 integer :: basis_size0,dimekb0,iekb,ii,ij,il,ilm,ilmn,iln,iln0,im,ios,iproj,iproj0,itypat,itypat0
 integer :: jekb,jlmn,jln,lmnmax0,lmn_size0,lmn2_size0,lpsang,nekb,nproj,npsang,pspso0,shape_type0
 integer :: usepaw0,vrspsp8
 real(dp) :: rpaw0,rshape0
 character(len=12) :: string
 character(len=500) :: message
!arrays
 integer,allocatable :: i1(:),i2(:),nprj(:),orbitals(:)
 real(dp),allocatable :: dij0(:),ekb0(:,:)

! *********************************************************************

!Check psddb8 version number (vrsio8) against DDB version number (vrsddb)
 if (vrsio8/=vrsddb) then
   write(message, '(a,i10,a,a,i10,a)' )&
&   'the psddb8 DDB version number=',vrsio8,ch10,&
&   'is not equal to the calling code DDB version number=',vrsddb,'.'
   MSG_WARNING(message)
 end if

!Check the value of choice
 if (choice<=0.or.choice>=3) then
   write(message, '(a,a,a,i10,a)' )&
&   'The permitted values for choice are 1 or 2.',ch10,&
&   'The calling routine asks ',choice,'.'
   MSG_BUG(message)
 end if

!==================================================================================
!First option: read psp characteristic from file =================================
!==================================================================================
 if (choice==1) then

   read(nunit,*)
   read(nunit, '(a12)' )string
   fullinit=1 ; if (lmnmax>0) indlmn(:,:,:)=0

!  --------------------------------------------
!  -----  NEW FORMAT (NCPP+PAW) ---------------
!  --------------------------------------------
   if (string=='  Descriptio')then
!    ==============================
!    ======== Common data =========
!    ==============================
     read(nunit, '(32x,i6)' )vrspsp8
     if (vrspsp8==vrsio8_old.or.vrspsp8==vrsio8_old_old) then
       usepaw0=0
       read(nunit, '(10x,i3,14x,i3,11x,i3)', iostat=ios)dimekb0,lmnmax0,usepaw0
       if(ios/=0)then
         backspace(nunit)
         read (nunit, '(10x,i3,14x,i3)' )dimekb0,lmnmax0
         usepaw0=0
       end if
     else if (vrspsp8==vrsio8) then
       read(nunit, '(10x,i3)') usepaw0
       if (usepaw/=usepaw0) then
         write(message, '(a,i1,a,i1,a)' )&
&         'usepaw is announced to be ',usepaw,' but read usepaw is ',usepaw0,' !'
         MSG_ERROR(message)
       end if
       if (usepaw==0) then
         read (nunit, '(10x,i3,14x,i3)' )dimekb0,lmnmax0
       end if
     end if

!    ==============================
!    === Norm-conserving psps =====
!    ==============================
     if (usepaw==0) then
       ekb(:,:)=zero
       ABI_MALLOC(ekb0,(dimekb,dimekb))
       do itypat=1,ntypat
         read(nunit, '(13x,i4,9x,i3,8x,i4)' )itypat0,pspso0,nekb
!        Check the compatibility with the main code dimensioning
         if(nekb>dimekb)then
           write(message, '(a,i8,a,a,a,i3,a)' )&
&           '  ',nekb,' components of ekb are announced',ch10,&
&           'but dimekb=',dimekb,'.'
           MSG_BUG(message)
         end if
         read(nunit,*)
         ilmn=0;iproj0=0
         do iekb=1,nekb
           read(nunit, '(3i6,3x,8d15.7)' ) iln,lpsang,iproj,(ekb0(ii,iekb),ii=1,min(nekb,4))
           if(nekb>4)then
             do jekb=5,nekb,4
               read(nunit, '(21x,8d15.7)' )(ekb0(ii,iekb),ii=jekb,min(nekb,jekb+3))
             end do
           end if
           if (lpsang==0.and.iproj>iproj0) iproj0=iproj
           if (useylm==1) then
             do im=-lpsang,lpsang
               ilmn=ilmn+1
               indlmn(1,ilmn,itypat)=lpsang
               indlmn(2,ilmn,itypat)=im
               indlmn(3,ilmn,itypat)=iproj
               indlmn(4,ilmn,itypat)=lpsang**2+lpsang+1+im
               indlmn(5,ilmn,itypat)=iln
               indlmn(6,ilmn,itypat)=1
               if (pspso0/=1.and.iln>(nekb-iproj0)/2) indlmn(6,ilmn,itypat)=2
             end do
           else
             ilmn=ilmn+1
             indlmn(1,ilmn,itypat)=lpsang
             indlmn(2,ilmn,itypat)=lpsang
             indlmn(3,ilmn,itypat)=iproj
             indlmn(4,ilmn,itypat)=lpsang**2+lpsang+1
             indlmn(5,ilmn,itypat)=iln
             indlmn(6,ilmn,itypat)=1
             if (pspso0/=1.and.iln>(nekb-iproj0)/2) indlmn(6,ilmn,itypat)=2
           end if
!          For the time being, only diagonal ekb are treated in abinit v3
           ekb(iekb,itypat)=ekb0(iekb,iekb)
!          For non-diagonal ekb, one could use:
!          do jekb=iekb to nekb
!          ekb(jekb+iekb*(iekb-1)/2,itypat)=ekb0(jekb,iekb)
!          end do
         end do
       end do
       ABI_FREE(ekb0)

!      ==============================
!      ============ PAW =============
!      ==============================
     else
       do itypat=1,ntypat
         read(nunit, '(12x,i4,12x,i3,12x,i5)' )itypat0,basis_size0,lmn_size0
         lmn2_size0=lmn_size0*(lmn_size0+1)/2
         ABI_MALLOC(orbitals,(basis_size0))
         read(nunit, '(20x,50i2)' ) orbitals(1:basis_size0)
         read(nunit, '(11x,f6.3,13x,i2,11x,f6.3)' ) rpaw0,shape_type0,rshape0
         read(nunit,'(24x,i3)') nekb
         read(nunit,*)
         ABI_MALLOC(dij0,(nekb))
         ABI_MALLOC(i1,(nekb))
         ABI_MALLOC(i2,(nekb))
         do ii=1,nekb,4
           read(nunit,'(3x,4(1x,i4,1x,i4,1x,d12.5))') (i1(ij),i2(ij),dij0(ij),ij=ii,min(ii+3,nekb))
         end do
         if (lmn_size0>lmnmax) then
           write(message, '(a,i5,3a,i5,a)' )&
&           'max. value of ',lmnmax,' for lmn_size is announced',ch10,&
&           'but ',lmn_size0,' is read.'
           MSG_BUG(message)
         end if
         if (allocated(pawtab(itypat)%dij0)) then
           if (lmn_size0>pawtab(itypat)%lmn_size) then
             write(message, '(a,i5,3a,i5,a)' )&
&             'lmn_size=,',pawtab(itypat)%lmn_size,' is announced',ch10,&
&             'but ',lmn_size0,' is read.'
             MSG_BUG(message)
           end if
         end if
         ABI_MALLOC(nprj,(0:maxval(orbitals)))
         ilmn=0;nprj=0
         do iln=1,basis_size0
           il=orbitals(iln)
           nprj(il)=nprj(il)+1
           do ilm=1,2*il+1
             indlmn(1,ilmn+ilm,itypat)=il
             indlmn(2,ilmn+ilm,itypat)=ilm-(il+1)
             indlmn(3,ilmn+ilm,itypat)=nprj(il)
             indlmn(4,ilmn+ilm,itypat)=il*il+ilm
             indlmn(5,ilmn+ilm,itypat)=iln
             indlmn(6,ilmn+ilm,itypat)=1
           end do
           ilmn=ilmn+2*il+1
         end do
         pawtab(itypat)%basis_size=basis_size0
         pawtab(itypat)%lmn_size  =lmn_size0
         pawtab(itypat)%lmn2_size =lmn2_size0
         pawtab(itypat)%shape_type=shape_type0
         pawtab(itypat)%rpaw      =rpaw0
         pawtab(itypat)%rshp      =rshape0
         if (.not.allocated(pawtab(itypat)%dij0))  then
           ABI_MALLOC(pawtab(itypat)%dij0,(lmn2_size0))
         end if
         pawtab(itypat)%dij0(1:lmn2_size0)=zero
         do ii=1,nekb
           ij=i1(ii)+i2(ii)*(i2(ii)-1)/2
           pawtab(itypat)%dij0(ij)=dij0(ii)
         end do
         ABI_FREE(nprj)
         ABI_FREE(orbitals)
         ABI_FREE(dij0)
         ABI_FREE(i1)
         ABI_FREE(i2)
       end do

     end if ! NCPP or PAW

!    --------------------------------------------
!    -----  OLD FORMAT (NCPP only) --------------
!    --------------------------------------------
   else if (string==' Description')then
     if (usepaw==1) then
       MSG_BUG("old DDB pspformat not compatible with PAW")
     end if

     read (nunit, '(10x,i3,10x,i3)' )nproj,npsang
     nekb=nproj*npsang
!    Check the compatibility with the main code dimensioning
     if(nekb>dimekb)then
       write(message, '(a,i8,a,a,a,i3,a)' )&
&       '  ',nekb,' components of ekb are announced',ch10,&
&       'but the maximum is dimekb=',dimekb,'.'
       MSG_BUG(message)
     end if
     if(useylm/=0)then
       MSG_BUG('useylm must be 0 !')
     end if
!    Read the data
     ABI_MALLOC(ekb0,(dimekb,dimekb))
     ekb0(:,:)=zero
     do itypat=1,ntypat
       read (nunit, '(13x,i4)' )ij
       do iproj=1,nproj
         read (nunit, '(6x,3d22.14)' )&
&         (ekb0(iproj+nproj*(ii-1),iproj+nproj*(ii-1)),ii=1,min(npsang,3))
         if(npsang>3)read (nunit, '(6x,3d22.14)' )&
&         (ekb0(iproj+nproj*(ii-1),iproj+nproj*(ii-1)),ii=4,npsang)
         do ii=1,npsang
           iekb=iproj+nproj*(ii-1)
           indlmn(1,iekb,itypat)=ii-1
           indlmn(2,iekb,itypat)=ii-1
           indlmn(3,iekb,itypat)=iproj
           indlmn(4,iekb,itypat)=ii**2-ii+1
           indlmn(5,iekb,itypat)=iekb
           indlmn(6,iekb,itypat)=1
!          For the time being, only diagonal ekb are treated in abinit v3
           ekb(iekb,itypat)=ekb0(iekb,iekb)
         end do
       end do
     end do
     ABI_FREE(ekb0)

!    --------------------------------------------
!    -----  OTHER CASES -------------------------
!    --------------------------------------------
   else if(string==' No informat')then
     fullinit=0
   else
     MSG_BUG('Error when reading the psp information')
   end if

!  Now, the number of blocks
   read(nunit,*)
   read(nunit,*)
   read(nunit, '(24x,i4)' )nblok

!  ==================================================================================
!  Second option: read psp characteristic from file ================================
!  ==================================================================================
 else if(choice==2)then

   write(nunit, '(a)' )' '
   if (fullinit==0)then
!    This possibility is used when the DDB is initialized,
!    and the ekb s are not available from the GS input file...
     write(nunit, '(a)' )' No information on the potentials yet '
   else

!    ==============================
!    === Norm-conserving psps =====
!    ==============================
     if (usepaw==0) then
       write(nunit, '(a)' )'  Description of the potentials (KB energies)'
       write(nunit, '(a,i6)' )'  vrsio8 (for pseudopotentials)=',vrsio8
       write(nunit, '(a,i3)' ) '  usepaw =',usepaw
       write(nunit, '(a,i3,a,i3,a,i3)' )'  dimekb =',dimekb,'       lmnmax=',lmnmax
       ABI_MALLOC(ekb0,(dimekb,dimekb))
       do itypat=1,ntypat
!        Compute nekb
         nekb=0
         do jlmn=1,lmnmax
           jln=indlmn(5,jlmn,itypat)
           if(jln>nekb)then
             nekb=jln
           end if
         end do
         write(nunit, '(a,i4,a,i3,a,i4)' ) &
&         '  Atom type= ',itypat,'   pspso=',pspso(itypat),'   nekb=',nekb
         write(nunit, '(a)' ) '  iln lpsang iproj  ekb(:)'
         iln0=0
         ekb0(:,:)=zero
         do ilmn=1,lmnmax
           iln =indlmn(5,ilmn,itypat)
           if (iln>iln0) then
             iln0=iln
             lpsang=indlmn(1,ilmn,itypat)
             iproj=indlmn(3,ilmn,itypat)
!            For the time being, only diagonal ekb are treated in abinit v3
             ekb0(iln,iln)=ekb(iln,itypat)
!            For non-diagonal ekb, one could use:
!            do ii=iln to nekb
!            ekb0(ii,iln)=ekb(ii+iln*(iln-1)/2,itypat)
!            end do
             write(nunit, '(3i6,3x,4es15.7)' ) iln,lpsang,iproj,(ekb0(ii,iln),ii=1,min(nekb,4))
             if(nekb>4)then
               do iekb=5,nekb,4
                 write(nunit, '(21x,4es15.7)' )(ekb0(ii,iekb),ii=iekb,min(nekb,iekb+3))
               end do
             end if
           end if
         end do
       end do
       ABI_FREE(ekb0)

!      ==============================
!      ============ PAW =============
!      ==============================
     else
       write(nunit, '(a)' )'  Description of the PAW dataset(s)'
       write(nunit, '(a,i6)' )'  vrsio8 (for pseudopotentials)=',vrsio8
       write(nunit, '(a,i3)' ) '  usepaw =',usepaw
       do itypat=1,ntypat
         iln0=0
         ABI_MALLOC(orbitals,(pawtab(itypat)%basis_size))
         do ilmn=1,pawtab(itypat)%lmn_size
           iln =indlmn(5,ilmn,itypat)
           if (iln>iln0) then
             iln0=iln;orbitals(iln)=indlmn(1,ilmn,itypat)
           end if
         end do
         write(nunit, '(a,i4,a,i3,a,i5)' ) &
&         '  Atom type=',itypat,' basis_size=',pawtab(itypat)%basis_size,&
&         '   lmn_size=',pawtab(itypat)%lmn_size
         write(nunit, '(a,50i2)' ) &
&         '    Basis functions=',orbitals(1:pawtab(itypat)%basis_size)
         write(nunit, '(a,f6.3,a,i2,a,f6.3)' ) &
&         '    r_PAW= ',pawtab(itypat)%rpaw,' shape_type= ',pawtab(itypat)%shape_type,&
&         '  r_shape= ',pawtab(itypat)%rshp
         nekb=0
         ABI_MALLOC(dij0,(pawtab(itypat)%lmn2_size))
         ABI_MALLOC(i1,(pawtab(itypat)%lmn2_size))
         ABI_MALLOC(i2,(pawtab(itypat)%lmn2_size))
         do jlmn=1,pawtab(itypat)%lmn_size
           ij=jlmn*(jlmn-1)/2
           do ilmn=1,jlmn
             if (abs(pawtab(itypat)%dij0(ij+ilmn))>tol16) then
               nekb=nekb+1;i1(nekb)=ilmn;i2(nekb)=jlmn
               dij0(nekb)=pawtab(itypat)%dij0(ij+ilmn)
             end if
           end do
         end do
         write(nunit,'(a,i3,a)') '    Dij0=     (only the ',nekb,' values different from zero)'
         write(nunit,'(2a)') '       i    j     Dij0        i    j     Dij0 ',&
&         '       i    j     Dij0        i    j     Dij0'
         do ii=1,nekb,4
           write(nunit,'(3x,4(1x,i4,1x,i4,1x,es12.5))') (i1(ij),i2(ij),dij0(ij),ij=ii,min(ii+3,nekb))
         end do
         ABI_FREE(dij0)
         ABI_FREE(i1)
         ABI_FREE(i2)
         ABI_FREE(orbitals)
       end do

     end if ! NCPP or PAW
   end if ! fullinit==0

!  Now, write the number of blocks
   write(nunit, '(a)' )' '
   write(nunit, '(a)' )' **** Database of total energy derivatives ****'
   write(nunit, '(a,i4)' ) ' Number of data blocks= ',nblok
 end if

end subroutine psddb8
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/ioddb8_in
!!
!! NAME
!! ioddb8_in
!!
!! FUNCTION
!! Open Derivative DataBase, and read preliminary information.
!! Note: only one processor reads the DDB.
!!
!! INPUTS
!! character(len=*)= name of input file
!! matom=maximum number of atoms
!! mband=maximum number of bands
!! mkpt=maximum number of special points
!! msym=maximum number of symetries
!! mtypat=maximum number of atom types
!! unddb=unit number for input
!! vrsddb=6 digit integer giving date, in form yymmdd for month=mm(1-12),
!!  day=dd(1-31), and year=yy(90-99 for 1990 to 1999,00-89 for 2000 to 2089),
!!  of current DDB version.
!!
!! OUTPUT
!! acell(3)=length scales of primitive translations (bohr)
!! amu(mtypat)=mass of the atoms (atomic mass unit)
!! dilatmx=the maximal dilatation factor
!! ecut=kinetic energy planewave cutoff (hartree)
!! ecutsm=smearing energy for plane wave kinetic energy (Ha)
!! intxc=control xc quadrature
!! iscf=parameter controlling scf or non-scf choice
!! ixc=exchange-correlation choice parameter
!! kpt(3,mkpt)=k point set (reduced coordinates)
!! kptnrm=normalisation of k points
!! natom=number of atoms in the unit cell
!! nband(mkpt)=number of bands at each k point, for each polarization
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!! nkpt=number of k points
!! nspden=number of spin-density components
!! nspinor=number of spinorial components of the wavefunctions
!! nsppol=1 for unpolarized, 2 for spin-polarized
!! nsym=number of symmetry elements in space group
!! ntypat=number of atom types
!! occ(mband*mkpt)=occupation number for each band and k
!! occopt=option for occupancies
!! pawecutdg=cut-off for fine "double grid" used in PAW calculations (unused for NCPP)
!! rprim(3,3)=dimensionless primitive translations in real space
!! dfpt_sciss=scissor shift (Ha)
!! spinat(3,matom)=initial spin of each atom, in unit of hbar/2
!! symafm(msym)=(anti)ferromagnetic part of symmetry operations
!! symrel(3,3,msym)=symmetry operations in real space
!! tnons(3,msym)=nonsymmorphic translations for symmetry operations
!! tolwfr=tolerance on largest wf residual
!! tphysel="physical" electronic temperature with FD occupations
!! tsmear=smearing width (or temperature) in Hartree
!! typat(matom)=type of each atom
!! usepaw=flag for PAW
!! wtk(mkpt)=weight assigned to each k point
!! xred(3,matom)=reduced atomic coordinates
!! zion(mtypat)=valence charge of each type of atom
!! znucl(mtypat)=atomic number of atom type
!!
!! PARENTS
!!      m_ddb,mblktyp1,mblktyp5,thmeig
!!
!! CHILDREN
!!
!! SOURCE

subroutine ioddb8_in(filnam,matom,mband,mkpt,msym,mtypat,unddb,vrsddb,&
&  acell,amu,dilatmx,ecut,ecutsm,intxc,iscf,ixc,kpt,kptnrm,&
&  natom,nband,ngfft,nkpt,nspden,nspinor,nsppol,nsym,ntypat,occ,occopt,&
&  pawecutdg,rprim,dfpt_sciss,spinat,symafm,symrel,tnons,tolwfr,tphysel,tsmear,&
&  typat,usepaw,wtk,xred,zion,znucl)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ioddb8_in'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: matom,mband,mkpt,msym,mtypat,unddb,vrsddb,usepaw
 integer,intent(out) :: intxc,iscf,ixc,natom,nkpt,nspden,nspinor,nsppol,nsym,ntypat,occopt
 real(dp),intent(out) :: dilatmx,ecut,ecutsm,pawecutdg,kptnrm,dfpt_sciss,tolwfr,tphysel,tsmear
 character(len=*),intent(in) :: filnam
!arrays
 integer,intent(out) :: nband(mkpt),ngfft(18),symafm(msym),symrel(3,3,msym),typat(matom)
 real(dp),intent(out) :: acell(3),amu(mtypat),kpt(3,mkpt),occ(mband*mkpt)
 real(dp),intent(out) :: rprim(3,3),spinat(3,matom),tnons(3,msym),wtk(mkpt)
 real(dp),intent(out) :: xred(3,matom),zion(mtypat),znucl(mtypat)

!Local variables -------------------------
!Set routine version number here:
!scalars
 integer,parameter :: vrsio8=100401,vrsio8_old=010929,vrsio8_old_old=990527
 integer :: bantot,ddbvrs,iband,ii,ij,ikpt,iline,im,usepaw0
 logical :: ddbvrs_is_current_or_old,testn,testv
 character(len=500) :: message
 character(len=6) :: name_old
!arrays
 character(len=12) :: name(9)

! *********************************************************************

!Check ioddb8 version number (vrsio8) against mkddb version number (vrsddb)
 if (vrsio8/=vrsddb) then
   write(message, '(a,i10,a,a,i10,a)' )&
&   'The input/output DDB version number=',vrsio8,ch10,&
&   'is not equal to the DDB version number=',vrsddb,'.'
   MSG_WARNING(message)
 end if

!Open the input derivative database.
 write(message,'(a,a)')' About to open file ',TRIM(filnam)
 call wrtout(std_out,message,'COLL')
 if (open_file(filnam,message,unit=unddb,form="formatted",status="old",action="read") /= 0) then
   MSG_ERROR(message)
 end if

!Check the compatibility of the input DDB with the DDB code
 read (unddb,*)
 read (unddb,*)
 read (unddb, '(20x,i10)' )ddbvrs

 !write(std_out,'(a,i10)')' ddbvrs=',ddbvrs
 if(ddbvrs/=vrsio8 .and. ddbvrs/=vrsio8_old .and. ddbvrs/=vrsio8_old_old)then
   write(message, '(a,i10,2a,3(a,i10),a)' )&
&   'The input DDB version number=',ddbvrs,' does not agree',ch10,&
&   'with the allowed code DDB version numbers,',vrsio8,', ',vrsio8_old,' and ',vrsio8_old_old,' .'
   MSG_BUG(message)
 end if

!Read the 4 n-integers, also testing the names of data, and checking that their value is acceptable.
!This is important to insure that any array has a sufficient dimension.
 read (unddb,*)
 read (unddb,*)
 read (unddb,*)
 testn=.true.; testv=.true.
 ddbvrs_is_current_or_old=(ddbvrs==vrsio8.or.ddbvrs==vrsio8_old)

!1. usepaw
 if(ddbvrs==vrsio8)then
   read (unddb, '(1x,a9,i10)' )name(1),usepaw0
 else
   usepaw0=0;name(1)='   usepaw'
 end if
 if(name(1)/='   usepaw')testn=.false.
 if(usepaw0/=usepaw)testv=.false.
!2. natom
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,i10)' )name(2),natom
 else
   read (unddb, '(1x,a6,i10)' )name_old,natom ; name(2)='   '//name_old
 end if
 if(name(2)/='    natom')testn=.false.
 if(natom<=0.or.natom>matom)testv=.false.
!3. nkpt
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,i10)' )name(3),nkpt
 else
   read (unddb, '(1x,a6,i10)' )name_old,nkpt ; name(3)='   '//name_old
 end if
 if(name(3)/='     nkpt')testn=.false.
 if(nkpt <=0.or.nkpt >mkpt )testv=.false.
!4. nsppol
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,i10)' )name(4),nsppol
 else
   read (unddb, '(1x,a6,i10)' )name_old,nsppol ; name(4)='   '//name_old
 end if
 if(name(4)/='   nsppol')testn=.false.
 if(nsppol<=0.or.nsppol>2)testv=.false.
!5. nsym
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,i10)' )name(5),nsym
 else
   read (unddb, '(1x,a6,i10)' )name_old,nsym ; name(5)='   '//name_old
 end if
 if(name(5)/='     nsym')testn=.false.
 if(nsym <=0.or.nsym >msym )testv=.false.
!6. ntypat
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,i10)' )name(6),ntypat
 else
   read (unddb, '(1x,a6,i10)' )name_old,ntypat ; name(6)='   '//name_old
 end if
 if(name(6)/='   ntypat' .and. name(6)/='    ntype')testn=.false.
 if(ntypat<=0.or.ntypat>mtypat)testv=.false.
!7. occopt
!Before reading nband, the last parameters that define
!the dimension of some array, need to know what is their
!representation, given by occopt
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,i10)' )name(7),occopt
 else
   read (unddb, '(1x,a6,i10)' )name_old,occopt ; name(7)='   '//name_old
 end if
 if(name(7)/='   occopt')testn=.false.
 if(occopt<0.or.occopt>8)testv=.false.
!Message if the names or values are not right
 if (.not.testn.or..not.testv) then
   write(message, '(a,a,a)' )' ioddb8_in : An error has been found in one',ch10,&
&   ' of the positive n-integers contained in the DDB : '
   call wrtout(std_out,message,'COLL')
   write(message, '(a)' )&
&   '               Expected                      Found     '
   call wrtout(std_out,message,'COLL')
   write(message, '(a,i10,a,a,a,i10)' )&
&   '    usepaw equal to   ',usepaw,'    ',trim(name(1)),' =',usepaw0
   call wrtout(std_out,message,'COLL')
   write(message, '(a,i10,a,a,a,i10)' )&
&   '    natom , lower than',matom+1,'    ',trim(name(2)),' =',natom
   call wrtout(std_out,message,'COLL')
   write(message, '(a,i10,a,a,a,i10)' )&
&   '    nkpt  , lower than',mkpt+1 ,'    ',trim(name(3)),' =',nkpt
   call wrtout(std_out,message,'COLL')
   write(message, '(a,i10,a,a,a,i10)' )&
&   '    nsppol, lower than',3      ,'    ',trim(name(4)),' =',nsppol
   call wrtout(std_out,message,'COLL')
   write(message, '(a,i10,a,a,a,i10)' )&
&   '    nsym  , lower than',msym+1 ,'    ',trim(name(5)),' =',nsym
   call wrtout(std_out,message,'COLL')
   write(message, '(a,i10,a,a,a,i10)' )&
&   '    ntypat, lower than',mtypat+1,'   ',trim(name(6)),' =',ntypat
   call wrtout(std_out,message,'COLL')
   write(message, '(a,a,a,i10)' )&
&   '    occopt,  between 0 and 7        ',trim(name(7)),' =',occopt
   call wrtout(std_out,message,'COLL')

   MSG_ERROR('See the error message above.')
 end if

!One more set of parameters define the dimensions of the
!array : nband. Morever, it depends on occopt and nkpt, and has to be
!tested after the test on nkpt is performed.
!8. nband
 if(occopt==2)then
   im=12
   do iline=1,(nkpt+11)/12
     name(:) = "" ! reset all line strings
     if(iline==(nkpt+11)/12)im=nkpt-12*(iline-1)
     if(ddbvrs_is_current_or_old)then
       read (unddb, '(1x,a9,5x,12i5)' )name(1),(nband((iline-1)*12+ii),ii=1,im)
     else
       read (unddb, '(1x,a6,5x,12i5)' )name_old,(nband((iline-1)*12+ii),ii=1,im) ; name(1)='   '//name_old
     end if
     if (iline==1) then
       call ddb_chkname(name(1),'    nband')
     else
       call ddb_chkname(name(1),'         ')
     end if
   end do
 else
   name(:) = "" ! reset all line strings
   if(ddbvrs_is_current_or_old)then
     read (unddb, '(1x,a9,i10)' )name(1),nband(1)
   else
     read (unddb, '(1x,a6,i10)' )name_old,nband(1) ; name(1)='   '//name_old
   end if
   call ddb_chkname(name(1),'    nband')
   if(nkpt>1)then
     do ikpt=2,nkpt
       nband(ikpt)=nband(1)
     end do
   end if
 end if

!check all nband values, and sum them
 bantot=0
 do ikpt=1,nkpt
   if(nband(ikpt)<0)then
     write(message, '(a,i4,a,i4,3a)' )&
&     'For ikpt = ',ikpt,'  nband = ',nband(ikpt),' is negative.',ch10,&
&     'Action: correct your DDB.'
     MSG_ERROR(message)
   else if(nband(ikpt)>mband)then
     write(message, '(a,i4,a,i4,a,a,i4,3a)' )&
&     'For ikpt = ',ikpt,', nband = ',nband(ikpt),ch10,&
&     'is larger than mband = ',mband,'.',ch10,&
&     'Action: recompile the calling code with a larger mband.'
     MSG_ERROR(message)
   end if
   bantot=bantot+nband(ikpt)
 end do

!Read the rest of variables, with check of the names
!9. acell
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,3d22.14)' )name(1),acell
 else
   read (unddb, '(1x,a6,3d22.14)' )name_old,acell ; name(1)='   '//name_old
 end if
 call ddb_chkname(name(1),'    acell')
!9. amu
 im=3
 do iline=1,(ntypat+2)/3
   if(iline==(ntypat+2)/3)im=ntypat-3*(iline-1)
   if(ddbvrs_is_current_or_old)then
     read (unddb, '(1x,a9,3d22.14)' )name(1),(amu((iline-1)*3+ii),ii=1,im)
   else
     read (unddb, '(1x,a6,3d22.14)' )name_old,(amu((iline-1)*3+ii),ii=1,im) ; name(1)='   '//name_old
   end if
   if (iline==1) then
     call ddb_chkname(name(1),'      amu')
   else
     call ddb_chkname(name(1),'         ')
   end if
 end do
!11. dilatmx
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,d22.14)' )name(1),dilatmx
   call ddb_chkname(name(1),'  dilatmx')
 else
   dilatmx=one
 end if
!12. ecut
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,d22.14)' )name(1),ecut
 else
   read (unddb, '(1x,a6,d22.14)' )name_old,ecut ; name(1)='   '//name_old
 end if
 call ddb_chkname(name(1),'     ecut')
!12b. pawecutdg (PAW only)
 if(ddbvrs==vrsio8.and.usepaw==1) then
   read (unddb, '(1x,a9,d22.14)' )name(1),pawecutdg
 else
   pawecutdg=ecut;name(1)='pawecutdg'
 end if
 call ddb_chkname(name(1),'pawecutdg')
!13. ecutsm
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,d22.14)' )name(1),ecutsm
   call ddb_chkname(name(1),'   ecutsm')
 else
   ecutsm=zero
 end if
!14. intxc
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,i10)' )name(1),intxc
   call ddb_chkname(name(1),'    intxc')
 else
   intxc=1
 end if
!15. iscf
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,i10)' )name(1),iscf
 else
   read (unddb, '(1x,a6,i10)' )name_old,iscf ; name(1)='   '//name_old
 end if
 call ddb_chkname(name(1),'     iscf')
!16. ixc
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,i10)' )name(1),ixc
 else
   read (unddb, '(1x,a6,i10)' )name_old,ixc ; name(1)='   '//name_old
 end if
 call ddb_chkname(name(1),'      ixc')
!17. kpt
 do iline=1,nkpt
   if(ddbvrs_is_current_or_old)then
     read (unddb, '(1x,a9,3d22.14)' )name(1),(kpt(ii,iline),ii=1,3)
   else
     read (unddb, '(1x,a6,3d22.14)' )name_old,(kpt(ii,iline),ii=1,3) ; name(1)='   '//name_old
   end if
   if (iline==1) then
     call ddb_chkname(name(1),'      kpt')
   else
     call ddb_chkname(name(1),'         ')
   end if
 end do
!18. kptnrm
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,d22.14)' )name(1),kptnrm
 else
   read (unddb, '(1x,a6,d22.14)' )name_old,kptnrm ; name(1)='   '//name_old
 end if
 call ddb_chkname(name(1),'   kptnrm')
!19. ngfft
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,5x,3i5)' )name(1),ngfft(1:3)
 else
   read (unddb, '(1x,a6,5x,3i5)' )name_old,ngfft(1:3) ; name(1)='   '//name_old
 end if
!For the time being, do not check the validity of the name,
!in order to accept both ng and ngfft
!20. nspden
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,i10)' )name(1),nspden
   call ddb_chkname(name(1),'   nspden')
 else
   nspden=0
 end if
!21. nspinor
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,i10)' )name(1),nspinor
   call ddb_chkname(name(1),'  nspinor')
 else
   nspinor=0
 end if
!22. occ
 if(occopt==2)then
   im=3
   do iline=1,(bantot+2)/3
     if(iline==(bantot+2)/3)im=bantot-3*(iline-1)
     if(ddbvrs_is_current_or_old)then
       read (unddb, '(1x,a9,3d22.14)' )name(1),(occ((iline-1)*3+ii),ii=1,im)
     else
       read (unddb, '(1x,a6,3d22.14)' )name_old,(occ((iline-1)*3+ii),ii=1,im) ; name(1)='   '//name_old
     end if
     if (iline==1) then
       call ddb_chkname(name(1),'      occ')
     else
       call ddb_chkname(name(1),'         ')
     end if
   end do
 else
   im=3
   do iline=1,(nband(1)+2)/3
     if(iline==(nband(1)+2)/3)im=nband(1)-3*(iline-1)
     if(ddbvrs_is_current_or_old)then
       read (unddb, '(1x,a9,3d22.14)' )name(1),(occ((iline-1)*3+ii),ii=1,im)
     else
       read (unddb, '(1x,a6,3d22.14)' )name_old,(occ((iline-1)*3+ii),ii=1,im) ; name(1)='   '//name_old
     end if
     if (iline==1) then
       call ddb_chkname(name(1),'      occ')
     else
       call ddb_chkname(name(1),'         ')
     end if
   end do
   if(nkpt>1)then
     do ikpt=2,nkpt
       do iband=1,nband(1)
         occ(iband+nband(1)*(ikpt-1))=occ(iband)
       end do
     end do
   end if
 end if
!23. rprim
 do iline=1,3
   if(ddbvrs_is_current_or_old)then
     read (unddb, '(1x,a9,3d22.14)' )name(1),(rprim(ii,iline),ii=1,3)
   else
     read (unddb, '(1x,a6,3d22.14)' )name_old,(rprim(ii,iline),ii=1,3) ; name(1)='   '//name_old
   end if
   if (iline==1) then
     call ddb_chkname(name(1),'    rprim')
   else
     call ddb_chkname(name(1),'         ')
   end if
 end do
!24. dfpt_sciss
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a11,d22.14)' )name(1),dfpt_sciss
   call ddb_chkname(name(1),'dfpt_sciss', 'sciss')
 else
   read (unddb, '(1x,a6,d22.14)' )name_old,dfpt_sciss ; name(1)=name_old
   call ddb_chkname(name(1),'sciss')
 end if
!25. spinat
 if(ddbvrs_is_current_or_old)then
   do iline=1,natom
     read (unddb, '(1x,a9,3d22.14)' )name(1),(spinat(ii,iline),ii=1,3)
     if (iline==1) then
       call ddb_chkname(name(1),'   spinat')
     else
       call ddb_chkname(name(1),'         ')
     end if
   end do
 else
!  spinat is set to zero by default in mrgddb.f
!  spinat(:,1:natom)=zero
 end if
!26. symafm
 if(ddbvrs_is_current_or_old)then
   im=12
   do iline=1,(nsym+11)/12
     if(iline==(nsym+11)/12)im=nsym-12*(iline-1)
     read (unddb, '(1x,a9,5x,12i5)' )name(1),(symafm((iline-1)*12+ii),ii=1,im)
     if (iline==1) then
       call ddb_chkname(name(1),'   symafm')
     else
       call ddb_chkname(name(1),'         ')
     end if
   end do
 else
!  symafm is set to 1 by default in mrgddb.f
!  symafm(1:nsym)=1
 end if
!27. symrel
 do iline=1,nsym
   if(ddbvrs_is_current_or_old)then
     read (unddb, '(1x,a9,5x,9i5)' )name(1),((symrel(ii,ij,iline),ii=1,3),ij=1,3)
   else
     read (unddb, '(1x,a6,5x,9i5)' )name_old,&
&     ((symrel(ii,ij,iline),ii=1,3),ij=1,3) ; name(1)='   '//name_old
   end if
   if (iline==1) then
     call ddb_chkname(name(1),'   symrel')
   else
     call ddb_chkname(name(1),'         ')
   end if
 end do
!28old. xred
 if(.not.ddbvrs_is_current_or_old)then
   do iline=1,natom
     read (unddb, '(1x,a6,3d22.14)' )name(1),(xred(ii,iline),ii=1,3)
   end do
!  No check of name, to allow the old tn
 end if
!28. tnons
 do iline=1,nsym
   if(ddbvrs_is_current_or_old)then
     read (unddb, '(1x,a9,3d22.14)' )name(1),(tnons(ii,iline),ii=1,3)
   else
     read (unddb, '(1x,a6,3d22.14)' )name_old,(tnons(ii,iline),ii=1,3) ; name(1)='   '//name_old
   end if
   if (iline==1) then
     call ddb_chkname(name(1),'    tnons')
   else
     call ddb_chkname(name(1),'         ')
   end if
 end do
!29. tolwfr
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,d22.14)' )name(1),tolwfr
 end if
!Do not check the name, in order to allow both tolwfr and wftol
!30. tphysel
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,d22.14)' )name(1),tphysel
   call ddb_chkname(name(1),'  tphysel')
 else
   tphysel=zero
 end if
!31. tsmear
 if(ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a9,d22.14)' )name(1),tsmear
   call ddb_chkname(name(1),'   tsmear')
 else
   tsmear=zero
 end if
!32. typat
 im=12
 do iline=1,(natom+11)/12
   if(iline==(natom+11)/12)im=natom-12*(iline-1)
   if(ddbvrs_is_current_or_old)then
     read (unddb, '(1x,a9,5x,12i5)' )name(1),(typat((iline-1)*12+ii),ii=1,im)
   else
     read (unddb, '(1x,a6,5x,12i5)' )name_old,(typat((iline-1)*12+ii),ii=1,im) ; name(1)='   '//name_old
   end if
   if (iline==1) then
!    Both type and typat are allowed => no check
!    call ddb_chkname(name(1),'    typat')
   else
     call ddb_chkname(name(1),'         ')
   end if
 end do
!33old. tolwfr
 if(.not.ddbvrs_is_current_or_old)then
   read (unddb, '(1x,a6,d22.14)' )name(1),tolwfr
 end if
!Do not check the name, in order to allow both tolwfr and wftol
!33. wtk
 im=3
 do iline=1,(nkpt+2)/3
   if(iline==(nkpt+2)/3)im=nkpt-3*(iline-1)
   if(ddbvrs_is_current_or_old)then
     read (unddb, '(1x,a9,3d22.14)' )name(1),(wtk((iline-1)*3+ii),ii=1,im)
   else
     read (unddb, '(1x,a6,3d22.14)' )name_old,(wtk((iline-1)*3+ii),ii=1,im) ; name(1)='   '//name_old
   end if
   if (iline==1) then
     call ddb_chkname(name(1),'      wtk')
   else
     call ddb_chkname(name(1),'         ')
   end if
 end do
!34. xred
 if(ddbvrs_is_current_or_old)then
   do iline=1,natom
     read (unddb, '(1x,a9,3d22.14)' )name(1),(xred(ii,iline),ii=1,3)
     if (iline==1) then
       call ddb_chkname(name(1),'     xred')
     else
       call ddb_chkname(name(1),'         ')
     end if
   end do
 end if
!35. znucl
 if(ddbvrs_is_current_or_old)then
   im=3
   do iline=1,(ntypat+2)/3
     if(iline==(ntypat+2)/3)im=ntypat-3*(iline-1)
     read (unddb, '(1x,a9,3d22.14)' )name(1),(znucl((iline-1)*3+ii),ii=1,im)
     if (iline==1) then
       call ddb_chkname(name(1),'    znucl')
     else
       call ddb_chkname(name(1),'         ')
     end if
   end do
 else
!  znucl is set to zero by default in mrgddb.f
!  znucl(:)=zero
 end if
!36. zion
 im=3
 do iline=1,(ntypat+2)/3
   if(iline==(ntypat+2)/3)im=ntypat-3*(iline-1)
   if(ddbvrs_is_current_or_old)then
     read (unddb, '(1x,a9,3d22.14)' )name(1),(zion((iline-1)*3+ii),ii=1,im)
   else
     read (unddb, '(1x,a6,3d22.14)' )name_old,(zion((iline-1)*3+ii),ii=1,im) ; name(1)='   '//name_old
   end if
   if (iline==1) then
!    Do not check the names, to allow both zion and znucl - the latter for 990527 format
!    call ddb_chkname(name(1),'     zion')
   else
     call ddb_chkname(name(1),'         ')
   end if
 end do

end subroutine ioddb8_in
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/rdddb9
!!
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
!! occopt=occupation option
!! rmet(3,3)=metric tensor in real space (bohr^2)
!! rprim(3,3)= primitive translation vectors
!! symrec(3,3,nsym)=3x3 matrices of the group symmetries (reciprocal space)
!! symrel(3,3,nsym)=3x3 matrices of the group symmetries (real space)
!! symafm(nsym)=Anti-ferromagnetic symmetries.
!! tnons(3,nsym)=fractional nonsymmorphic translations
!! typat(natom)=type integer for each atom in cell
!! ucvol=unit cell volume in bohr**3
!! usepaw= 0 for non paw calculation; =1 for paw calculation
!! xcart(3,natom)=atomic cartesian coordinates
!! xred(3,natom)=fractional dimensionless atomic coordinates
!! zion(ntypat)=charge on each type of atom (real number)
!! znucl(ntypat)=Nuclear charge for each type of pseudopotential
!!
!! PARENTS
!!      m_ddb
!!
!! CHILDREN
!!
!! SOURCE

subroutine rdddb9(acell,atifc,amu,ddb,&
& ddbun,dimekb,filnam,gmet,gprim,indsym,iout,&
& lmnmax,mband,mpert,msize,msym,&
& natifc,natom,nkpt,nsym,ntypat,&
& occopt,rmet,rprim,symrec,symrel,symafm,&
& tnons,typat,ucvol,usepaw,xcart,xred,zion,znucl)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rdddb9'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_41_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
! NOTE: these are used for dimensioning and then re-assigned in ioddb8.
!   This is almost definitely bad practice. In particular
!    it should be indsym(4,msym,natom),
!   and
!    the allocation allocate(kpt(3,nkpt)) is strange
!scalars
 integer,intent(in) :: ddbun,dimekb,iout,lmnmax,mband,mpert,msize,msym,natifc
 integer,intent(inout) :: natom,nkpt,nsym,ntypat,occopt,usepaw
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
 integer :: mtypat,mkpt,matom
 integer :: choice,fullinit,iblok,intxc,iscf,isym,ixc
 integer :: nsize,nspden,nspinor,nsppol,nunit,timrev,useylm,vrsddb
 real(dp),parameter :: tolsym8=tol8
 real(dp) :: dilatmx,ecut,ecutsm,kptnrm,pawecutdg,dfpt_sciss,tolwfr
 real(dp) :: tphysel,tsmear
 character(len=500) :: message
!arrays
 integer :: ngfft(18),symq(4,2,msym)
 integer,allocatable :: car3flg(:,:,:,:,:,:),carflg(:,:,:,:),indlmn(:,:,:)
 integer,allocatable :: nband(:),pspso(:),tmpflg(:,:,:,:,:,:)
 real(dp) :: gprimd(3,3),qpt(3),rprimd(3,3)
 real(dp),allocatable :: d2cart(:,:,:,:,:),d3cart(:,:,:,:,:,:,:),ekb(:,:)
 real(dp),allocatable :: kpt(:,:),occ(:),spinat(:,:),tmpval(:,:,:,:,:,:,:)
 real(dp),allocatable :: wtk(:)
 type(pawtab_type),allocatable :: pawtab(:)

! *********************************************************************

 DBG_ENTER("COLL")

!Read the DDB information
 vrsddb=DDB_VERSION

!The checking of pseudopotentials is not done presently so that dimensions are fake
 ABI_MALLOC(ekb,(dimekb,ntypat))
 ABI_MALLOC(indlmn,(6,lmnmax,ntypat))
 ABI_MALLOC(pspso,(ntypat))
 ABI_DATATYPE_ALLOCATE(pawtab,(ntypat*usepaw))
 call pawtab_nullify(pawtab)

 ABI_MALLOC(kpt,(3,nkpt))
 ABI_MALLOC(nband,(nkpt))
 ABI_MALLOC(occ,(nkpt*mband*msppol))
 ABI_MALLOC(spinat,(3,natom))
 ABI_MALLOC(wtk,(nkpt))

!Open the input derivative database file and read the preliminary information
!Note that in this call, mkpt has been replaced by nkpt, mtypat by ntypat, and matom by natom.
 nunit=ddbun

! To avoid aliasing
 matom = natom
 mkpt = nkpt
 mtypat = ntypat

 call ioddb8_in(filnam,matom,mband,mkpt,msym,mtypat,nunit,vrsddb,&
& acell,amu,dilatmx,ecut,ecutsm,intxc,iscf,ixc,kpt,kptnrm,&
& natom,nband,ngfft,nkpt,nspden,nspinor,nsppol,nsym,ntypat,occ,occopt,&
& pawecutdg,rprim,dfpt_sciss,spinat,symafm,symrel,tnons,tolwfr,tphysel,tsmear,&
& typat,usepaw,wtk,xred,zion,znucl)

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

!Read the psp information of the input DDB
 useylm=usepaw;choice=1
 call psddb8 (choice,dimekb,ekb,fullinit,indlmn,lmnmax,&
& ddb%nblok,ntypat,nunit,pawtab,pspso,usepaw,useylm,vrsddb)

!Check the correctness of some input parameters, and perform small treatment if needed.
 call chkin9(atifc,natifc,natom)

!Read the blocks from the input database, and close it.
 write(message, '(3a,i0,a)' )ch10,ch10,' rdddb9: read ',ddb%nblok,' blocks from the input DDB '
 call wrtout(std_out,message,'COLL')
 nunit=ddbun

 do iblok=1,ddb%nblok
   call read_blok8(ddb,iblok,mband,mpert,msize,nkpt,nunit)

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
     call d2sym3(tmpflg,tmpval,indsym,mpert,natom,nsym,qpt,symq,symrec,symrel,timrev)

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

     nsize=3*mpert*3*mpert*3*mpert
     ABI_MALLOC(tmpflg,(3,mpert,3,mpert,3,mpert))
     ABI_MALLOC(tmpval,(2,3,mpert,3,mpert,3,mpert))

     tmpflg(:,:,:,:,:,:) = reshape(ddb%flg(1:nsize,iblok), shape = (/3,mpert,3,mpert,3,mpert/))
     tmpval(1,:,:,:,:,:,:) = reshape(ddb%val(1,1:nsize,iblok), shape = (/3,mpert,3,mpert,3,mpert/))
     tmpval(2,:,:,:,:,:,:) = reshape(ddb%val(2,1:nsize,iblok), shape = (/3,mpert,3,mpert,3,mpert/))

     call d3sym(tmpflg,tmpval,indsym,mpert,natom,nsym,symrec,symrel)

     ABI_MALLOC(d3cart,(2,3,mpert,3,mpert,3,mpert))
     ABI_MALLOC(car3flg,(3,mpert,3,mpert,3,mpert))

     call nlopt(tmpflg,car3flg,tmpval,d3cart,gprimd,mpert,natom,rprimd,ucvol)

     ddb%flg(1:nsize,iblok) = reshape(car3flg, shape = (/3*mpert*3*mpert*3*mpert/))
     ddb%val(1,1:nsize,iblok) = reshape(d3cart(1,:,:,:,:,:,:), shape = (/3*mpert*3*mpert*3*mpert/))
     ddb%val(2,1:nsize,iblok) = reshape(d3cart(2,:,:,:,:,:,:), shape = (/3*mpert*3*mpert*3*mpert/))

     ABI_FREE(d3cart)
     ABI_FREE(car3flg)
     ABI_FREE(tmpflg)
     ABI_FREE(tmpval)
   end if
 end do ! iblok

 close(ddbun)

 write(message,'(a)' )' Now the whole DDB is in central memory '
 call wrtout(std_out,message,'COLL')
 call wrtout(iout,message,'COLL')

 ABI_FREE(ekb)
 ABI_FREE(indlmn)
 ABI_FREE(kpt)
 ABI_FREE(nband)
 ABI_FREE(occ)
 ABI_FREE(pspso)
 ABI_FREE(spinat)
 ABI_FREE(wtk)

 call pawtab_free(pawtab)
 ABI_DATATYPE_DEALLOCATE(pawtab)

 DBG_EXIT("COLL")

end subroutine rdddb9
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/chkin9
!!
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
!!
!! SOURCE

subroutine chkin9(atifc,natifc,natom)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chkin9'
!End of the abilint section

 implicit none

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
   write(message, '(a,i0,a,a,a,i0,a,a,a)' )&
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
       write(message, '(a,i0,a,a,a,a,a,i0,a,a,a)' )&
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
!!  gprimd(3,3)=dimensional primitive translations for
!!              reciprocal space(bohr^-1)
!!  mpert =maximum number of ipert
!!  natom= number of atoms
!!  rprimd(3,3)=dimensional primitive translations (bohr)
!!  ucvol=unit cell volume (bohr^3)
!!
!! OUTPUT
!! carflg(3,mpert,3,mpert,3,mpert)=1 if the element of d3cart has been
!!   calculated, 0 otherwise
!! d3cart(2,3,mpert,3,mpert,3,mpert)=matrix of third-order energy
!!   derivatives in cartesian coordinates
!!
!! PARENTS
!!      m_ddb,nonlinear
!!
!! CHILDREN
!!
!! SOURCE

subroutine nlopt(blkflg,carflg,d3,d3cart,gprimd,mpert,natom,rprimd,ucvol)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nlopt'
!End of the abilint section

 implicit none

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
!!  brav = 1 -> simple lattice; 2 -> face-centered cubic;
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
!!      anaddb,dfpt_looppert,eph,m_effective_potential_file
!!
!! CHILDREN
!!
!! SOURCE

subroutine ddb_from_file(ddb,filename,brav,natom,natifc,atifc,crystal,comm,prtvol)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ddb_from_file'
 use interfaces_41_geometry
!End of the abilint section

 implicit none

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
!arrays
 integer,allocatable :: symrec(:,:,:),symrel(:,:,:),symafm(:),indsym(:,:,:),typat(:)
 real(dp) :: acell(3),gmet(3,3),gprim(3,3),rmet(3,3),rprim(3,3),rprimd(3,3)
 real(dp),allocatable :: amu(:),xcart(:),xred(:,:),zion(:),znucl(:),tnons(:,:)
 character(len=132),allocatable :: title(:)
 character(len=500) :: message

! ************************************************************************

 DBG_ENTER("COLL")

! Must read natom from the DDB before being able to allocate some arrays needed for invars9
 call ddb_getdims(dimekb,filename,lmnmax,mband,mtyp,msym,ddb_natom,nblok,nkpt,ntypat,get_unit(),usepaw,DDB_VERSION,comm)
 if (ddb_natom /= natom) then
   MSG_ERROR(sjoin("input natom:",itoa(natom),"does not agree with DDB value:",itoa(natom)))
 end if

 mpert=natom+6
 msize=3*mpert*3*mpert; if (mtyp==3) msize=msize*3*mpert

 ! Allocate arrays depending on msym
 ! (which is actually fixed to nsym inside inprep8)
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

   call ddb_malloc(ddb,msize,nblok,natom,ntypat)

   call rdddb9(acell,atifc,amu,ddb,&
&   ddbun,dimekb,filename,gmet,gprim,indsym,ab_out,&
&   lmnmax,mband,mpert,msize,msym,&
&   natifc,ddb_natom,nkpt,nsym,ntypat,&
&   occopt,rmet,rprim,symrec,symrel,symafm,&
&   tnons,typat,ucvol,usepaw,xcart,xred,zion,znucl)

   close(ddbun)

   ABI_FREE(symrec)
   ABI_FREE(indsym)
   ABI_FREE(xcart)

   ! Renormalize rprim to possibly satisfy the constraint abs(rprim(1,2))=half when brav/=1
   ! This section is needed to preserver the behaviour of the old implementation.
   if (brav/=1 .and. abs(abs(rprim(1,2))-half)>tol10) then
     if(abs(rprim(1,2))<tol6)then
       write(message, '(a,i0,7a)' )&
&       'The input DDB value of brav is ',brav,',',ch10,&
&       'and the one of rprim(1,2) is zero.',ch10,&
&       'These are incompatible',ch10,&
&       'Action: check the value of brav and rprim(1,2) in your DDB.'
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
   ! The information on it is contained in the four arrays
   !   ddb%flg(msize,nblok) : blok flag for each element
   !   ddb%qpt(9,nblok)     : blok wavevector (unnormalized)
   !   ddb%nrm(3,nblok)     : blok wavevector normalization
   !   ddb%typ(nblok)       : blok type
 end if

 if (xmpi_comm_size(comm) > 1) then
   call ddb_bcast (ddb, master, comm)
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

!Initialize crystal_t object.
 call mkrdim(acell,rprim,rprimd)

!FIXME: These variables are hardcoded
 npsp = ntypat; space_group = 0; timrev = 2
 use_antiferro=.FALSE. !;  use_antiferro=(nspden==2.and.nsppol==1)
 ABI_MALLOC(title, (ntypat))

 do ii=1,ntypat
   write(title(ii),'(a,i0)')"No title for typat ",ii
 end do

!Warning znucl is dimension with ntypat = nspsp hence alchemy is not supported here
 call crystal_init(ddb%amu,Crystal,space_group,natom,npsp,ntypat,nsym,rprimd,typat,xred,&
&  zion,znucl,timrev,use_antiferro,.FALSE.,title,&
&  symrel=symrel,tnons=tnons,symafm=symafm)

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
!!
!! SOURCE

subroutine carttransf(blkflg,blkval2,carflg,gprimd,iqpt,mband,&
& mpert,msize,natom,nblok,nkpt,rprimd)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'carttransf'
!End of the abilint section

 implicit none

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
!!
!! SOURCE

subroutine carteig2d(blkflg,blkval,carflg,d2cart,&
& gprimd,iblok,mpert,natom,nblok,rprimd)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'carteig2d'
!End of the abilint section

 implicit none

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
!!
!! OUTPUT
!! zeff(3,3,natom)=effective charge on each atom, versus electric
!!  field and atomic displacement. Note the following convention:
!!  zeff(electric field direction, atomic direction, atom number)
!! dielt(3,3)=dielectric tensor
!!
!! PARENTS
!!      m_ddb
!!
!! CHILDREN
!!
!! SOURCE

subroutine dtech9(blkval,dielt,iblok,mpert,natom,nblok,zeff)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dtech9'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: iblok,mpert,natom,nblok
!arrays
 real(dp),intent(in) :: blkval(2,3,mpert,3,mpert,nblok)
 real(dp),intent(out) :: dielt(3,3),zeff(3,3,natom)

!Local variables -------------------------
!scalars
 integer :: depl,elec,elec1,elec2,iatom
 character(len=1000) :: message

! *********************************************************************

!Extration of effectives charges
 do iatom=1,natom
   do elec=1,3
     do depl=1,3
       zeff(elec,depl,iatom)=0.5*&
&       (blkval(1,depl,iatom,elec,natom+2,iblok)+&
&       blkval(1,elec,natom+2,depl,iatom,iblok))
     end do
   end do
 end do

!Extration of dielectric tensor
 do elec1=1,3
   do elec2=1,3
     dielt(elec1,elec2)=blkval(1,elec1,natom+2,elec2,natom+2,iblok)
   end do
 end do

 write(message,'(a,3es16.6,3es16.6,3es16.6)' )&
& ' Dielectric Tensor ',&
& dielt(1,1),dielt(1,2),dielt(1,3),&
& dielt(2,1),dielt(2,2),dielt(2,3),&
& dielt(3,1),dielt(3,2),dielt(3,3)


 call wrtout(std_out,message,'COLL')

 write(message,'(a)' ) ' Effectives Charges '
 call wrtout(std_out,message,'COLL')
 do iatom=1,natom
   write(message,'(a,i4,3es16.6,3es16.6,3es16.6)' )' atom ',iatom,&
&   zeff(1,1,iatom),zeff(1,2,iatom),zeff(1,3,iatom),&
&   zeff(2,1,iatom),zeff(2,2,iatom),zeff(2,3,iatom),&
&   zeff(3,1,iatom),zeff(3,2,iatom),zeff(3,3,iatom)
    call wrtout(std_out,message,'COLL')
 end do

end subroutine dtech9
!!***

!----------------------------------------------------------------------

!!****f* m_ddb/dtchi
!!
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
!!
!! SOURCE

subroutine dtchi(blkval,dchide,dchidt,mpert,natom,ramansr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dtchi'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: mpert,natom,ramansr
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
 do iatom = 1, natom
   do depl = 1,3
     do elfd1 = 1,3
       do elfd2 = 1,3
         dchidt(iatom,depl,elfd1,elfd2) = d3cart(1,depl,iatom,elfd1,natom+2,elfd2,natom+2)
       end do
     end do
   end do
 end do

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

!DEBUG
!sumrule(:,:,:) = 0._dp
!do elfd2 = 1,3
!do elfd1 = 1,3
!do depl = 1,3
!do iatom = 1, natom
!sumrule(depl,elfd1,elfd2) = sumrule(depl,elfd1,elfd2) + &
!&     dchidt(iatom,depl,elfd1,elfd2)
!end do
!end do
!end do
!end do
!do depl = 1,3
!write(ab_out,'(6x,i2,3(3x,f16.9))') depl,sumrule(depl,1,1:3)
!write(ab_out,'(8x,3(3x,f16.9))') sumrule(depl,2,1:3)
!write(ab_out,'(8x,3(3x,f16.9))') sumrule(depl,3,1:3)
!write(ab_out,*)
!end do
!ENDEBUG

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

integer function ddb_get_etotal(ddb,etotal) result(iblok)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ddb_get_etotal'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 real(dp),intent(out) :: etotal
 type(ddb_type),intent(in) :: ddb

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

 call gtblk9(ddb,iblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,rftyp)

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
!!    (0=> no selection, 1=> trace only, 2=> symmetric part only)                       !!
!!
!! SIDE EFFECTS
!!  ddb<type(ddb_type)>=
!!    The block with the effective charges is modified if charge neutrality is imposed.
!!
!! OUTPUT
!!  dielt(3,3) = Macroscopic dielectric tensor
!!  zeff(3,3,natom)=effective charge on each atom, versus electric field and atomic displacement
!!  iblok=Index of the block containing the data. 0 if block is not found.
!!
!! NOTES
!!  dielt and zeff are initialized to one_3D and zero if the derivatives are not available in the DDB file.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function ddb_get_dielt_zeff(ddb,crystal,rftyp,chneut,selectz,dielt,zeff) result(iblok)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ddb_get_dielt_zeff'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: rftyp,chneut,selectz
 type(ddb_type),intent(inout) :: ddb
 type(crystal_t),intent(in) :: crystal
!arrays
 real(dp),intent(out) :: dielt(3,3),zeff(3,3,crystal%natom)

!Local variables -------------------------
!scalars
 integer :: ii
 character(len=500) :: message
!arrays
 integer :: rfelfd(4),rfphon(4),rfstrs(4)
 real(dp) :: qphnrm(3),qphon(3,3)

! *********************************************************************

 ! Look for the Gamma Block in the DDB
 qphon(:,1)=zero
 qphnrm(1)=zero
 rfphon(1:2)=1
 rfelfd(1:2)=2
 rfstrs(1:2)=0

 !write(std_out,*)"ddb%mpert",ddb%mpert

 call gtblk9(ddb,iblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,rftyp)

 ! Compute effective charges and dielectric tensor only if the Gamma-blok was found in the DDB
 ! In case it was not found, iblok = 0
 zeff=zero; dielt=zero; dielt(1,1)=one; dielt(2,2)=one; dielt(3,3)=one

 if (iblok/=0) then
   write(message, '(a,a,(80a),a,a,a,a)' ) ch10,('=',ii=1,80),ch10,ch10,&
&   ' Dielectric Tensor and Effective Charges ',ch10
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

   ! Make the imaginary part of the Gamma block vanish
   write(message, '(a,a,a,a,a)'  ) ch10,&
&   ' anaddb : Zero the imaginary part of the Dynamical Matrix at Gamma,',ch10,&
&   '   and impose the ASR on the effective charges ',ch10
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

   ! Impose the charge neutrality on the effective charges and eventually select some parts of the effective charges
   call chneu9(chneut,ddb%val(:,:,iblok),ddb%mpert,ddb%natom,ddb%ntypat,selectz,Crystal%typat,Crystal%zion)

   ! Extraction of the dielectric tensor and the effective charges
   call dtech9(ddb%val,dielt,iblok,ddb%mpert,ddb%natom,ddb%nblok,zeff)
 end if ! iblok not found

end function ddb_get_dielt_zeff
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

integer function ddb_get_dchidet(ddb,ramansr,dchide,dchidt) result(iblok)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ddb_get_dchidet'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: ramansr
 type(ddb_type),intent(in) :: ddb
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
 rfphon(1)  = 1 ; rfphon(2:3) = 0
 rfelfd(:)  = 2
 rfstrs(:)  = 0
 rftyp = 3

 call gtblk9(ddb,iblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,rftyp)

 if (iblok /= 0) then
   call dtchi(ddb%val(:,:,iblok),dchide,dchidt,ddb%mpert,ddb%natom,ramansr)
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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ddb_get_asrq0'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: asr,rftyp
 type(ddb_type),intent(inout) :: ddb
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

 call gtblk9(ddb,asrq0%iblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,rftyp)
 ! this is to maintain the old behaviour in which the arrays where allocated and set to zero in anaddb.
 ABI_MALLOC(asrq0%d2asr, (2,3,ddb%natom,3,ddb%natom))
 asrq0%d2asr = zero

 ! TODO: Tests with asr = 3,4  [v5][t83] and [v5][t84]
 ! fail if I don't allocated these arrays because the code
 ! is accessing the data without checking if the correction has been computed....
 dims = 3*ddb%natom*(3*ddb%natom-1) / 2
 ABI_CALLOC(asrq0%uinvers, (dims, dims))
 ABI_CALLOC(asrq0%vtinvers,(dims, dims))
 ABI_CALLOC(asrq0%singular, (dims))

 if (asrq0%iblok == 0) return
 iblok = asrq0%iblok

 select case (asr)
 case (0)
   continue

 case (1,2)
   call asria_calc(asr,asrq0%d2asr,ddb%val(:,:,iblok),ddb%mpert,ddb%natom)

 case (3,4)
   ! Rotational invariance for 1D and 0D systems
   ! Compute uinvers, vtinvers and singular matrices.
   !dims = 3*ddb%natom*(3*ddb%natom-1) / 2
   !ABI_CALLOC(asrq0%uinvers, (dims, dims))
   !ABI_CALLOC(asrq0%vtinvers,(dims, dims))
   !ABI_CALLOC(asrq0%singular, (dims))

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
!!
!! SOURCE

subroutine ddb_diagoq(ddb, crystal, qpt, asrq0, symdynmat, rftyp, phfrq, displ_cart, &
                      out_eigvec,out_displ_red)   ! Optional [out]


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ddb_diagoq'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: rftyp,symdynmat
 type(ddb_type),intent(in) :: ddb
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

 call gtblk9(ddb,iblok,qphon_padded,qphnrm,rfphon,rfelfd,rfstrs,rftyp)
 if (iblok == 0) then
   MSG_ERROR(sjoin("Cannot find q-point ", ktoa(qpt)," in DDB file"))
 end if

 ! Copy the dynamical matrix in d2cart
 d2cart(:,1:ddb%msize) = ddb%val(:,:,iblok)

 ! Eventually impose the acoustic sum rule based on previously calculated d2asr
 call asrq0_apply(asrq0, natom, ddb%mpert, ddb%msize, crystal%xcart, d2cart)

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
!!
!! SOURCE

subroutine asrq0_apply(asrq0, natom, mpert, msize, xcart, d2cart)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'asrq0_apply'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom, msize, mpert
 type(asrq0_t),intent(inout) :: asrq0
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
!!      anaddb,m_effective_potential,m_effective_potential_file
!!
!! CHILDREN
!!
!! SOURCE

subroutine asrq0_free(asrq0)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'asrq0_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(asrq0_t),intent(inout) :: asrq0

! ************************************************************************

 ! real
 if (allocated(asrq0%d2asr)) then
   ABI_FREE(asrq0%d2asr)
 end if

 if (allocated(asrq0%singular)) then
   ABI_FREE(asrq0%singular)
 end if

 if (allocated(asrq0%uinvers)) then
   ABI_FREE(asrq0%uinvers)
 end if

 if (allocated(asrq0%vtinvers)) then
   ABI_FREE(asrq0%vtinvers)
 end if

end subroutine asrq0_free
!!***

!!****f* m_ddb/ddb_chkname
!! NAME
!! ddb_chkname
!!
!! FUNCTION
!! This small subroutine check the identity of its argument,
!! who are a6 names, and eventually send a message and stop
!! if they are found unequal
!!
!! INPUTS
!! nmfond= name which has to be checked
!! nmxpct= name expected for nmfond
!! nmxpct2= eventual second optional name (backward compatibility)
!!
!! OUTPUT
!!
!! TODO
!! Describe the inputs
!!
!! PARENTS
!!      m_ddb
!!
!! CHILDREN
!!
!! SOURCE

subroutine ddb_chkname(nmfond,nmxpct,nmxpct2)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ddb_chkname'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 character(len=*),intent(in) :: nmfond,nmxpct
 character(len=*),intent(in),optional :: nmxpct2

!Local variables-------------------------------
!scalars
 logical :: found
 character(len=500) :: nmfond_,nmxpct_,nmxpct2_
 character(len=500) :: message

! *********************************************************************

 nmxpct_ = trim(adjustl(nmxpct))
 nmfond_ = trim(adjustl(nmfond))

 found = (nmxpct_ == nmfond_)

 if (present(nmxpct2) .and. .not. found) then
   nmxpct2_ = trim(adjustl(nmxpct2))
   found = (nmxpct2_==nmfond_)
 end if

 if (.not. found) then
   write(message, '(a,a,a,a,a,a,a,a,a,a,a)' )&
&   '  Reading DDB, expected name was "',trim(nmxpct_),'"',ch10,&
&   '               and name found is "',trim(nmfond_),'"',ch10,&
&   '  Likely your DDB is incorrect.',ch10,&
&   '  Action : correct your DDB, or contact the ABINIT group.'
   MSG_ERROR(message)
 end if

end subroutine ddb_chkname
!!***

!----------------------------------------------------------------------

END MODULE m_ddb

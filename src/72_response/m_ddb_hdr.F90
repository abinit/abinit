!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_ddb_hdr
!! NAME
!!  m_ddb_hdr
!!
!! FUNCTION
!!  This module contains the declaration of data types and methods
!!  to handle the header of the DDB files.
!!
!! COPYRIGHT
!! Copyright (C) 2011-2019 ABINIT group (MJV, XG, MT, MM, MVeithen, MG, PB, JCC, GA)
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

MODULE m_ddb_hdr

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors
 use m_abicore
 use m_xmpi
 use m_dtset

 use m_copy,      only : alloc_copy
 use m_pawtab,    only : pawtab_type, pawtab_nullify, pawtab_free !, pawtab_copy
 use m_psps,      only : psps_copy, psps_free
 use m_io_tools,  only : open_file
 use m_copy,      only : alloc_copy
 use m_fstrings,  only : sjoin

 implicit none

 private

 public :: ddb_getdims      ! Open a DDB file and read basic dimensions and variables.
 public :: ioddb8_in        ! Temporary
 public :: psddb8           ! Temporary

 type,public :: ddb_hdr_type

   integer :: ddb_version   ! Version of the DDB file

   integer :: matom
   integer :: mband
   integer :: mkpt
   integer :: msym
   integer :: mtypat

   integer :: intxc
   integer :: iscf
   integer :: ixc
   integer :: natom
   integer :: nkpt
   integer :: nspden
   integer :: nspinor
   integer :: nsppol
   integer :: nsym
   integer :: ntypat
   integer :: occopt
   integer :: usepaw

   integer :: nblok         ! Number of blocks in the ddb
   integer :: mblktyp       ! Max block type
                            ! 0 = Total energy
                            ! 1 = 2nd derivatives (non-stat.)
                            ! 2 = 2nd derivatives (stationary)
                            ! 3 = 3rd derivatives
                            ! 4 = 1st derivatives
                            ! 5 = 2nd eigenvalue derivatives
                            !
   integer :: fullinit      ! Whether the full info on the pseudo is present
                            ! TODO rename this variable

   real(dp) :: dilatmx
   real(dp) :: ecut
   real(dp) :: ecutsm
   real(dp) :: kptnrm
   real(dp) :: pawecutdg
   real(dp) :: dfpt_sciss
   real(dp) :: tolwfr
   real(dp) :: tphysel
   real(dp) :: tsmear

   character(len=fnlen) :: dscrpt

   integer :: ngfft(18)
   real(dp) :: acell(3)
   real(dp) :: rprim(3,3)

   integer,allocatable :: nband(:)
   ! nband(mkpt*nsppol)

   integer,allocatable :: symafm(:)
   ! symafm(msym)

   integer,allocatable :: symrel(:,:,:)
   ! symrel(3,3,msym)

   integer,allocatable :: typat(:)
   ! typat(matom)

   real(dp),allocatable :: amu(:)
   ! amu(mtypat)

   real(dp),allocatable :: kpt(:,:)
   ! kpt(3,mkpt)

   real(dp),allocatable :: occ(:)
   ! occ(mband*mkpt*nsppol)

   real(dp),allocatable :: spinat(:,:)
   ! spinat(3,matom)

   real(dp),allocatable :: tnons(:,:)
   ! tnons(3,msym)

   real(dp),allocatable :: wtk(:)
   ! wtk(mkpt)

   real(dp),allocatable :: xred(:,:)
   ! xred(3,matom)

   real(dp),allocatable :: zion(:)
   ! zion(mtypat)

   real(dp),allocatable :: znucl(:)
   ! znucl(mtypat)

   type(pawtab_type),allocatable :: pawtab(:)
   ! pawtab(psps%ntypat*psps%usepaw)

   type(pseudopotential_type) :: psps


 end type ddb_hdr_type

 public :: ddb_hdr_init            ! Construct object
 public :: ddb_hdr_malloc          ! Allocate dynamic memory.
 public :: ddb_hdr_free            ! Free dynamic memory.
 public :: ddb_hdr_open_write      ! Open the DDB file and write the header.
 public :: ddb_hdr_open_read       ! Open the DDB file and read the header.
 public :: ddb_hdr_compare         ! Compare two DDB headers.

CONTAINS  !===========================================================
!!***

!----------------------------------------------------------------------

!!****f* m_ddb_hdr/ddb_hdr_init
!! NAME
!! ddb_hdr_init
!!
!! FUNCTION
!!   Initialize a ddb_hdr object from a dataset.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      dfpt_looppert,eig2tot,gstate,nonlinear,respfn
!!
!! CHILDREN
!!
!! SOURCE

subroutine ddb_hdr_init(ddb_hdr, dtset, psps, pawtab, ddb_version, dscrpt, &
                        nblok, xred, occ, ngfft)

!Arguments ------------------------------------
 type(ddb_hdr_type),intent(out) :: ddb_hdr
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)
 character(len=*),intent(in) :: dscrpt
 integer,intent(in) :: ddb_version
 integer,intent(in) :: nblok
 integer,intent(in),optional :: ngfft(18)
 real(dp),intent(in),optional :: xred(3,dtset%natom)
 real(dp),intent(in),optional :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)

!Local variables -------------------------
 integer :: ii, nn


! ************************************************************************

 ddb_hdr%nblok = nblok
 ddb_hdr%ddb_version = ddb_version
 ddb_hdr%dscrpt = trim(dscrpt)
 if (present(ngfft)) then
   ddb_hdr%ngfft = ngfft
 else
   ddb_hdr%ngfft = dtset%ngfft
 end if

 call psps_copy(psps, ddb_hdr%psps)

 ! Copy scalars from dtset
 ddb_hdr%matom = dtset%natom
 ddb_hdr%natom = dtset%natom
 ddb_hdr%mband = dtset%mband
 !ddb_hdr%mkpt = dtset%maxnkpt
 ddb_hdr%mkpt = dtset%nkpt
 ddb_hdr%nkpt = dtset%nkpt
 !ddb_hdr%msym = dtset%maxnsym
 ddb_hdr%msym = dtset%nsym
 ddb_hdr%nsym = dtset%nsym
 ddb_hdr%mtypat = dtset%ntypat
 ddb_hdr%ntypat = dtset%ntypat

 ddb_hdr%nspden = dtset%nspden
 ddb_hdr%nspinor = dtset%nspinor
 ddb_hdr%nsppol = dtset%nsppol

 ddb_hdr%occopt = dtset%occopt
 ddb_hdr%usepaw = dtset%usepaw

 ddb_hdr%intxc = dtset%intxc
 ddb_hdr%ixc = dtset%ixc
 ddb_hdr%iscf = dtset%iscf

 ddb_hdr%dilatmx = dtset%dilatmx
 ddb_hdr%ecut = dtset%ecut
 ddb_hdr%ecutsm = dtset%ecutsm
 ddb_hdr%pawecutdg = dtset%pawecutdg
 ddb_hdr%kptnrm = dtset%kptnrm
 ddb_hdr%dfpt_sciss = dtset%dfpt_sciss
 ddb_hdr%tolwfr = 1.0_dp  ! dummy
 ddb_hdr%tphysel = dtset%tphysel
 ddb_hdr%tsmear = dtset%tsmear

 ddb_hdr%fullinit = 1

 call ddb_hdr_malloc(ddb_hdr)

 ! Copy arrays from dtset
 ddb_hdr%acell(:) = dtset%acell_orig(1:3,1)
 ddb_hdr%rprim(:,:) = dtset%rprim_orig(1:3,1:3,1)
 ddb_hdr%amu(:) = dtset%amu_orig(:,1)
 ddb_hdr%nband(:) = dtset%nband(1:ddb_hdr%mkpt*ddb_hdr%nsppol)
 ddb_hdr%symafm(:) = dtset%symafm(1:ddb_hdr%msym)
 ddb_hdr%symrel(:,:,:) = dtset%symrel(1:3,1:3,1:ddb_hdr%msym)
 ddb_hdr%typat(:) = dtset%typat(1:ddb_hdr%matom)
 ddb_hdr%kpt(:,:) = dtset%kpt(1:3,1:ddb_hdr%mkpt)
 ddb_hdr%wtk(:) = dtset%wtk(1:ddb_hdr%mkpt)
 ddb_hdr%spinat(:,:) = dtset%spinat(1:3,1:ddb_hdr%matom)
 ddb_hdr%tnons(:,:) = dtset%tnons(1:3,1:ddb_hdr%msym)
 ddb_hdr%zion(:) = dtset%ziontypat(1:ddb_hdr%mtypat)
 ddb_hdr%znucl(:) = dtset%znucl(1:ddb_hdr%mtypat)

 if (present(xred)) then
   ddb_hdr%xred(:,:) = xred(1:3,1:ddb_hdr%matom)
 else
   ddb_hdr%xred(:,:) = dtset%xred_orig(1:3,1:ddb_hdr%matom,1)
 end if
 if (present(occ)) then
   ddb_hdr%occ(:) = occ(1:ddb_hdr%mband*ddb_hdr%mkpt*ddb_hdr%nsppol)
 else
   ddb_hdr%occ(:) = dtset%occ_orig(1:ddb_hdr%mband*ddb_hdr%mkpt*ddb_hdr%nsppol,1)
 end if


 ! GA: I had way too much problems implementing pawtab_copy.
 !     The script check-libpaw would report all sorts of errors.
 !     Therefore, I do a cheap copy here, copying only the relevant info.
 !call pawtab_copy(pawtab, ddb_hdr%pawtab)
 nn=size(pawtab)
 if (nn.gt.0) then
   do ii=1,nn
    ddb_hdr%pawtab(ii)%basis_size = pawtab(ii)%basis_size
    ddb_hdr%pawtab(ii)%lmn_size = pawtab(ii)%lmn_size
    ddb_hdr%pawtab(ii)%lmn2_size = pawtab(ii)%lmn2_size
    ddb_hdr%pawtab(ii)%rpaw = pawtab(ii)%rpaw
    ddb_hdr%pawtab(ii)%rshp = pawtab(ii)%rshp
    ddb_hdr%pawtab(ii)%shape_type = pawtab(ii)%shape_type
    if (allocated(pawtab(ii)%dij0)) then
      call alloc_copy(pawtab(ii)%dij0, ddb_hdr%pawtab(ii)%dij0)
    end if
   end do
 end if

end subroutine ddb_hdr_init
!!***

!!****f* m_ddb_hdr/ddb_hdr_malloc
!! NAME
!! ddb_hdr_malloc
!!
!! FUNCTION
!!  Allocate dynamic memory.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_ddb_hdr
!!
!! CHILDREN
!!
!! SOURCE

subroutine ddb_hdr_malloc(ddb_hdr)

!Arguments ------------------------------------
 type(ddb_hdr_type),intent(inout) :: ddb_hdr

! ************************************************************************

 ! integer
 ABI_MALLOC(ddb_hdr%nband,(ddb_hdr%mkpt*ddb_hdr%nsppol))
 ABI_MALLOC(ddb_hdr%symafm,(ddb_hdr%msym))
 ABI_MALLOC(ddb_hdr%symrel,(3,3,ddb_hdr%msym))
 ABI_MALLOC(ddb_hdr%typat,(ddb_hdr%matom))

 ! real
 ABI_MALLOC(ddb_hdr%amu,(ddb_hdr%mtypat))
 ABI_MALLOC(ddb_hdr%kpt,(3, ddb_hdr%mkpt))
 ABI_MALLOC(ddb_hdr%occ,(ddb_hdr%mband*ddb_hdr%mkpt*ddb_hdr%nsppol))
 ABI_MALLOC(ddb_hdr%spinat,(3, ddb_hdr%matom))
 ABI_MALLOC(ddb_hdr%tnons,(3, ddb_hdr%msym))
 ABI_MALLOC(ddb_hdr%wtk,(ddb_hdr%mkpt))
 ABI_MALLOC(ddb_hdr%xred,(3, ddb_hdr%matom))
 ABI_MALLOC(ddb_hdr%zion,(ddb_hdr%mtypat))
 ABI_MALLOC(ddb_hdr%znucl,(ddb_hdr%mtypat))

 ! types
 ABI_DATATYPE_ALLOCATE(ddb_hdr%pawtab,(ddb_hdr%psps%ntypat*ddb_hdr%psps%usepaw))
 call pawtab_nullify(ddb_hdr%pawtab)

end subroutine ddb_hdr_malloc
!!***

!----------------------------------------------------------------------

!!****f* m_ddb_hdr/ddb_hdr_free
!! NAME
!! ddb_hdr_free
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      anaddb,dfpt_looppert,eig2tot,gstate,m_ddb,m_effective_potential_file
!!      m_gruneisen,mblktyp1,mblktyp5,mrgddb,nonlinear,respfn,thmeig
!!
!! CHILDREN
!!
!! SOURCE

subroutine ddb_hdr_free(ddb_hdr)

!Arguments ------------------------------------
 type(ddb_hdr_type),intent(inout) :: ddb_hdr

! ************************************************************************

 ! integer
 ABI_SFREE(ddb_hdr%nband)
 ABI_SFREE(ddb_hdr%symafm)
 ABI_SFREE(ddb_hdr%symrel)
 ABI_SFREE(ddb_hdr%typat)

 ! real
 ABI_SFREE(ddb_hdr%amu)
 ABI_SFREE(ddb_hdr%kpt)
 ABI_SFREE(ddb_hdr%occ)
 ABI_SFREE(ddb_hdr%spinat)
 ABI_SFREE(ddb_hdr%tnons)
 ABI_SFREE(ddb_hdr%wtk)
 ABI_SFREE(ddb_hdr%xred)
 ABI_SFREE(ddb_hdr%zion)
 ABI_SFREE(ddb_hdr%znucl)

 ! types
 call psps_free(ddb_hdr%psps)

 if (allocated(ddb_hdr%pawtab)) then
   call pawtab_free(ddb_hdr%pawtab)
   ABI_DATATYPE_DEALLOCATE(ddb_hdr%pawtab)
 end if

end subroutine ddb_hdr_free
!!***

!----------------------------------------------------------------------

!!****f* m_ddb_hdr/ddb_hdr_open_write
!! NAME
!! ddb_hdr_open_write
!!
!! FUNCTION
!!  Open the DDB file and write the header.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      ddb_interpolate,dfpt_looppert,eig2tot,gstate,mblktyp1,mblktyp5
!!      nonlinear,respfn
!!
!! CHILDREN
!!
!! SOURCE

subroutine ddb_hdr_open_write(ddb_hdr, filnam, unddb, fullinit)

!Arguments ------------------------------------
 type(ddb_hdr_type),intent(inout) :: ddb_hdr
 character(len=*),intent(in) :: filnam
 integer,intent(in) :: unddb
 integer,intent(in),optional :: fullinit

!Local variables -------------------------
 integer :: choice

! ************************************************************************

 if (present(fullinit)) ddb_hdr%fullinit = fullinit

 call ddb_io_out(ddb_hdr%dscrpt,filnam,ddb_hdr%matom,ddb_hdr%mband,&
&  ddb_hdr%mkpt,ddb_hdr%msym,ddb_hdr%mtypat,unddb,ddb_hdr%ddb_version,&
&  ddb_hdr%acell,ddb_hdr%amu,ddb_hdr%dilatmx,ddb_hdr%ecut,ddb_hdr%ecutsm,&
&  ddb_hdr%intxc,ddb_hdr%iscf,ddb_hdr%ixc,ddb_hdr%kpt,ddb_hdr%kptnrm,&
&  ddb_hdr%natom,ddb_hdr%nband,ddb_hdr%ngfft,ddb_hdr%nkpt,ddb_hdr%nspden,&
&  ddb_hdr%nspinor,ddb_hdr%nsppol,ddb_hdr%nsym,ddb_hdr%ntypat,ddb_hdr%occ,&
&  ddb_hdr%occopt,ddb_hdr%pawecutdg,ddb_hdr%rprim,ddb_hdr%dfpt_sciss,&
&  ddb_hdr%spinat,ddb_hdr%symafm,ddb_hdr%symrel,ddb_hdr%tnons,ddb_hdr%tolwfr,&
&  ddb_hdr%tphysel,ddb_hdr%tsmear,ddb_hdr%typat,ddb_hdr%usepaw,ddb_hdr%wtk,&
&  ddb_hdr%xred,ddb_hdr%zion,ddb_hdr%znucl)


 choice=2  ! Write
 call psddb8(choice,ddb_hdr%psps%dimekb,ddb_hdr%psps%ekb,ddb_hdr%fullinit,&
&  ddb_hdr%psps%indlmn,ddb_hdr%psps%lmnmax,ddb_hdr%nblok,ddb_hdr%ntypat,unddb,&
&  ddb_hdr%pawtab,ddb_hdr%psps%pspso,ddb_hdr%psps%usepaw,ddb_hdr%psps%useylm,&
&  ddb_hdr%ddb_version)


end subroutine ddb_hdr_open_write
!!***

!----------------------------------------------------------------------

!!****f* m_ddb_hdr/ddb_hdr_open_read
!! NAME
!! ddb_hdr_open_read
!!
!! FUNCTION
!!  Open the DDB file and read the header.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      anaddb,m_ddb,m_effective_potential_file,m_gruneisen,mblktyp1,mblktyp5
!!      mrgddb,thmeig
!!
!! CHILDREN
!!
!! SOURCE

subroutine ddb_hdr_open_read(ddb_hdr, filnam, unddb, ddb_version, comm, &
&        matom,mtypat,mband,mkpt,msym,dimekb,lmnmax,usepaw,dimonly)

!Arguments ------------------------------------
 type(ddb_hdr_type),intent(inout) :: ddb_hdr
 character(len=*),intent(in) :: filnam
 integer,intent(in) :: unddb
 integer,intent(in) :: ddb_version
 integer,intent(in),optional :: comm
 integer,intent(in),optional :: matom,mtypat,mband,mkpt,msym
 integer,intent(in),optional :: dimekb,lmnmax,usepaw
 integer,intent(in),optional :: dimonly

!Local variables -------------------------
 integer :: choice
 integer :: mblktyp,nblok
 integer :: matom_l,mtypat_l,mband_l,mkpt_l,msym_l
 integer :: dimekb_l,lmnmax_l,usepaw_l
 integer :: comm_l

! ************************************************************************

 if (present(comm)) then
   comm_l = comm
 else
   comm_l = xmpi_comm_self
 end if

 ! Read the dimensions from the DDB
 call ddb_getdims(dimekb_l,filnam,lmnmax_l,mband_l,mblktyp, &
&       msym_l,matom_l,nblok,mkpt_l,mtypat_l,unddb,usepaw_l, &
&       ddb_version,comm_l)

 close(unddb)

 if (present(matom)) matom_l = matom
 if (present(mtypat)) mtypat_l = mtypat
 if (present(mband)) mband_l = mband
 if (present(mkpt)) mkpt_l = mkpt
 if (present(msym)) msym_l = msym
 if (present(dimekb)) dimekb_l = dimekb
 if (present(lmnmax)) lmnmax_l = lmnmax
 if (present(usepaw)) usepaw_l = usepaw

 ddb_hdr%ddb_version = ddb_version
 ddb_hdr%usepaw = usepaw_l
 ddb_hdr%mband = mband_l
 ddb_hdr%matom = matom_l
 ddb_hdr%msym = msym_l
 ddb_hdr%mtypat = mtypat_l
 ddb_hdr%mkpt = mkpt_l
 ddb_hdr%nsppol = 1     ! GA: Is nsppol not read?? Have to fix this...

 ddb_hdr%nblok = nblok
 ddb_hdr%mblktyp = mblktyp

 ddb_hdr%psps%dimekb = dimekb_l
 ddb_hdr%psps%ntypat = mtypat_l
 ddb_hdr%psps%lmnmax = lmnmax_l
 ddb_hdr%psps%usepaw = usepaw_l
 ddb_hdr%psps%useylm = usepaw_l  ! yep...

 ddb_hdr%natom = ddb_hdr%matom
 ddb_hdr%nkpt = ddb_hdr%mkpt
 ddb_hdr%nsym = ddb_hdr%msym
 ddb_hdr%ntypat = ddb_hdr%mtypat

 if (present(dimonly)) then
   if (dimonly==1) return
 end if

 ! Allocate the memory
 call ddb_hdr_malloc(ddb_hdr)

 ABI_ALLOCATE(ddb_hdr%psps%indlmn,(6,ddb_hdr%psps%lmnmax,ddb_hdr%mtypat))
 ABI_ALLOCATE(ddb_hdr%psps%pspso,(ddb_hdr%mtypat))
 ABI_ALLOCATE(ddb_hdr%psps%ekb,(ddb_hdr%psps%dimekb,ddb_hdr%mtypat))

 ! This is needed to read the DDBs in the old format
 ! GA : Not sure why this is required
 ddb_hdr%symafm(:)=1
 if(ddb_hdr%mtypat>=1)then
   ddb_hdr%psps%pspso(:)=0
   ddb_hdr%znucl(:)=zero
   ddb_hdr%psps%ekb(:,:)=zero
 end if
 if(ddb_hdr%matom>=1)then
   ddb_hdr%spinat(:,:)=zero
 end if


 ! Note: the maximum parameters (matom, mkpt, etc.) are inputs to ioddb8_in
 !       wile the actual parameters (natom, nkpt, etc.) are outputs
 call ioddb8_in(filnam,ddb_hdr%matom,ddb_hdr%mband,&
&       ddb_hdr%mkpt,ddb_hdr%msym,ddb_hdr%mtypat,unddb,ddb_version,&
&       ddb_hdr%acell,ddb_hdr%amu,ddb_hdr%dilatmx,ddb_hdr%ecut,ddb_hdr%ecutsm,&
&       ddb_hdr%intxc,ddb_hdr%iscf,ddb_hdr%ixc,ddb_hdr%kpt,ddb_hdr%kptnrm,&
&       ddb_hdr%natom,ddb_hdr%nband,ddb_hdr%ngfft,ddb_hdr%nkpt,ddb_hdr%nspden,&
&       ddb_hdr%nspinor,ddb_hdr%nsppol,ddb_hdr%nsym,ddb_hdr%ntypat,&
&       ddb_hdr%occ,ddb_hdr%occopt,ddb_hdr%pawecutdg,ddb_hdr%rprim,&
&       ddb_hdr%dfpt_sciss,ddb_hdr%spinat,ddb_hdr%symafm,ddb_hdr%symrel,&
&       ddb_hdr%tnons,ddb_hdr%tolwfr,ddb_hdr%tphysel,ddb_hdr%tsmear,&
&       ddb_hdr%typat,ddb_hdr%usepaw,ddb_hdr%wtk,ddb_hdr%xred,ddb_hdr%zion,&
&       ddb_hdr%znucl)


!  Read the psp information of the input DDB
   choice=1  ! Read
   call psddb8(choice,ddb_hdr%psps%dimekb,ddb_hdr%psps%ekb,ddb_hdr%fullinit,&
&       ddb_hdr%psps%indlmn,ddb_hdr%psps%lmnmax,&
&       ddb_hdr%nblok,ddb_hdr%ntypat,unddb,ddb_hdr%pawtab,ddb_hdr%psps%pspso,&
&       ddb_hdr%psps%usepaw,ddb_hdr%psps%useylm,ddb_version)

end subroutine ddb_hdr_open_read
!!***

!----------------------------------------------------------------------

!!****f* m_ddb_hdr/ddb_hdr_compare
!! NAME
!! ddb_hdr_compare
!!
!! FUNCTION
!!  Compare two DDB headers and raise error if they differ.
!!  Also, complete psps information if one has more info than the other.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      mblktyp1,mblktyp5
!!
!! CHILDREN
!!
!! SOURCE

subroutine ddb_hdr_compare(ddb_hdr1, ddb_hdr2)

!Arguments ------------------------------------
 type(ddb_hdr_type),intent(inout) :: ddb_hdr1, ddb_hdr2

!Local variables -------------------------

! ************************************************************************

!    Should also compare indlmn and pspso ... but suppose that
!    the checking of ekb is enough for the psps.
!    Should also compare many other variables ... this is still
!    to be done ...
 call compare_ddb_variables(&
&     ddb_hdr1%matom, ddb_hdr2%matom,&
&     ddb_hdr1%mtypat, ddb_hdr2%mtypat,&
&     ddb_hdr1%mkpt, ddb_hdr2%mkpt,&
&     ddb_hdr1%mband, ddb_hdr2%mband,&
&     ddb_hdr1%msym, ddb_hdr2%msym,&
&     ddb_hdr1%acell, ddb_hdr2%acell,&
&     ddb_hdr1%amu, ddb_hdr2%amu,&
&     ddb_hdr1%psps%dimekb,&
&     ddb_hdr1%ecut, ddb_hdr2%ecut,&
&     ddb_hdr1%psps%ekb, ddb_hdr2%psps%ekb,&
&     ddb_hdr1%fullinit, ddb_hdr2%fullinit,&
&     ddb_hdr1%iscf, ddb_hdr2%iscf,&
&     ddb_hdr1%ixc, ddb_hdr2%ixc,&
&     ddb_hdr1%kpt, ddb_hdr2%kpt,&
&     ddb_hdr1%kptnrm, ddb_hdr2%kptnrm,&
&     ddb_hdr1%natom, ddb_hdr2%natom,&
&     ddb_hdr1%nband, ddb_hdr2%nband,&
&     ddb_hdr1%ngfft, ddb_hdr2%ngfft,&
&     ddb_hdr1%nkpt, ddb_hdr2%nkpt,&
&     ddb_hdr1%nsppol, ddb_hdr2%nsppol,&
&     ddb_hdr1%nsym, ddb_hdr2%nsym,&
&     ddb_hdr1%ntypat, ddb_hdr2%ntypat,&
&     ddb_hdr1%occ, ddb_hdr2%occ,&
&     ddb_hdr1%occopt, ddb_hdr2%occopt,&
&     ddb_hdr1%pawecutdg, ddb_hdr2%pawecutdg,&
&     ddb_hdr1%pawtab, ddb_hdr2%pawtab,&
&     ddb_hdr1%rprim, ddb_hdr2%rprim,&
&     ddb_hdr1%dfpt_sciss, ddb_hdr2%dfpt_sciss,&
&     ddb_hdr1%symrel, ddb_hdr2%symrel,&
&     ddb_hdr1%tnons, ddb_hdr2%tnons,&
&     ddb_hdr1%tolwfr, ddb_hdr2%tolwfr,&
&     ddb_hdr1%typat, ddb_hdr2%typat,&
&     ddb_hdr1%usepaw,&
&     ddb_hdr1%wtk, ddb_hdr2%wtk,&
&     ddb_hdr1%xred, ddb_hdr2%xred,&
&     ddb_hdr1%zion, ddb_hdr2%zion)

end subroutine ddb_hdr_compare
!!***

!----------------------------------------------------------------------

!!****f* m_ddb_hdr/psddb8
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
!!      m_ddb_hdr
!!
!! CHILDREN
!!
!! SOURCE

subroutine psddb8 (choice,dimekb,ekb,fullinit,indlmn,lmnmax,&
&          nblok,ntypat,nunit,pawtab,pspso,usepaw,useylm,vrsddb)

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
    'the psddb8 DDB version number=',vrsio8,ch10,&
    'is not equal to the calling code DDB version number=',vrsddb,'.'
   MSG_WARNING(message)
 end if

!Check the value of choice
 if (choice<=0.or.choice>=3) then
   write(message, '(a,a,a,i10,a)' )&
    'The permitted values for choice are 1 or 2.',ch10,&
    'The calling routine asks ',choice,'.'
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
         write(message, '(a,i1,a,i1,a)' )'usepaw is announced to be ',usepaw,' but read usepaw is ',usepaw0,' !'
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
            '  ',nekb,' components of ekb are announced',ch10,'but dimekb=',dimekb,'.'
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
             'max. value of ',lmnmax,' for lmn_size is announced',ch10,'but ',lmn_size0,' is read.'
           MSG_BUG(message)
         end if
         if (allocated(pawtab(itypat)%dij0)) then
           if (lmn_size0>pawtab(itypat)%lmn_size) then
             write(message, '(a,i5,3a,i5,a)' )&
              'lmn_size=,',pawtab(itypat)%lmn_size,' is announced',ch10,'but ',lmn_size0,' is read.'
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
        '  ',nekb,' components of ekb are announced',ch10,'but the maximum is dimekb=',dimekb,'.'
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
         read (nunit, '(6x,3d22.14)' )(ekb0(iproj+nproj*(ii-1),iproj+nproj*(ii-1)),ii=1,min(npsang,3))
         if(npsang>3)read (nunit, '(6x,3d22.14)' )(ekb0(iproj+nproj*(ii-1),iproj+nproj*(ii-1)),ii=4,npsang)
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
         write(nunit, '(a,i4,a,i3,a,i4)' )'  Atom type= ',itypat,'   pspso=',pspso(itypat),'   nekb=',nekb
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
          '  Atom type=',itypat,' basis_size=',pawtab(itypat)%basis_size,&
          '   lmn_size=',pawtab(itypat)%lmn_size
         write(nunit, '(a,50i2)' ) &
          '    Basis functions=',orbitals(1:pawtab(itypat)%basis_size)
         write(nunit, '(a,f6.3,a,i2,a,f6.3)' ) &
          '    r_PAW= ',pawtab(itypat)%rpaw,' shape_type= ',pawtab(itypat)%shape_type,&
          '  r_shape= ',pawtab(itypat)%rshp
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
          '       i    j     Dij0        i    j     Dij0'
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

!!****f* m_ddb_hdr/ioddb8_in
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
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
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
!!      m_ddb_hdr
!!
!! CHILDREN
!!
!! SOURCE

subroutine ioddb8_in(filnam,matom,mband,mkpt,msym,mtypat,unddb,vrsddb,&
&  acell,amu,dilatmx,ecut,ecutsm,intxc,iscf,ixc,kpt,kptnrm,&
&  natom,nband,ngfft,nkpt,nspden,nspinor,nsppol,nsym,ntypat,occ,occopt,&
&  pawecutdg,rprim,dfpt_sciss,spinat,symafm,symrel,tnons,tolwfr,tphysel,tsmear,&
&  typat,usepaw,wtk,xred,zion,znucl)

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
    'The input/output DDB version number=',vrsio8,ch10,&
    'is not equal to the DDB version number=',vrsddb,'.'
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
    'The input DDB version number=',ddbvrs,' does not agree',ch10,&
    'with the allowed code DDB version numbers,',vrsio8,', ',vrsio8_old,' and ',vrsio8_old_old,' .'
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
    ' of the positive n-integers contained in the DDB : '
   call wrtout(std_out,message,'COLL')
   write(message, '(a)' )&
    '               Expected                      Found     '
   call wrtout(std_out,message,'COLL')
   write(message, '(a,i10,a,a,a,i10)' )&
    '    usepaw equal to   ',usepaw,'    ',trim(name(1)),' =',usepaw0
   call wrtout(std_out,message,'COLL')
   write(message, '(a,i10,a,a,a,i10)' )&
    '    natom , lower than',matom+1,'    ',trim(name(2)),' =',natom
   call wrtout(std_out,message,'COLL')
   write(message, '(a,i10,a,a,a,i10)' )&
    '    nkpt  , lower than',mkpt+1 ,'    ',trim(name(3)),' =',nkpt
   call wrtout(std_out,message,'COLL')
   write(message, '(a,i10,a,a,a,i10)' )&
    '    nsppol, lower than',3      ,'    ',trim(name(4)),' =',nsppol
   call wrtout(std_out,message,'COLL')
   write(message, '(a,i10,a,a,a,i10)' )&
    '    nsym  , lower than',msym+1 ,'    ',trim(name(5)),' =',nsym
   call wrtout(std_out,message,'COLL')
   write(message, '(a,i10,a,a,a,i10)' )&
    '    ntypat, lower than',mtypat+1,'   ',trim(name(6)),' =',ntypat
   call wrtout(std_out,message,'COLL')
   write(message, '(a,a,a,i10)' )&
    '    occopt,  between 0 and 7        ',trim(name(7)),' =',occopt
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

!!****f* m_ddb_hdr/ddb_getdims
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
!!      m_ddb_hdr
!!
!! CHILDREN
!!
!! SOURCE

subroutine ddb_getdims(dimekb,filnam,lmnmax,mband,mblktyp,msym,natom,nblok,nkpt,ntypat,unddb,usepaw,vrsddb,comm)

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

!!****f* m_ddb_hdr/inprep8
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
!!      m_ddb_hdr
!!
!! CHILDREN
!!
!! SOURCE


subroutine inprep8 (dimekb,filnam,lmnmax,mband,mblktyp,msym,natom,nblok,nkpt,&
& ntypat,unddb,usepaw,vrsddb)

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
     call wrtout(std_out,message,'COLL')  ! GA: why are we printing this?
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
   !usepaw=0   ! GA: usepaw is also declared earlier in the header
               !     and it is that earlier value that usepaw will
               !     be compared in ioddb8_in, so there is no reason
               !     to override the value here.

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

!!****f* m_ddb_hdr/ddb_chkname
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
!!      m_ddb_hdr
!!
!! CHILDREN
!!
!! SOURCE

subroutine ddb_chkname(nmfond,nmxpct,nmxpct2)

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
&   'Reading DDB, expected name was "',trim(nmxpct_),'"',ch10,&
&   '             and name found is "',trim(nmfond_),'"',ch10,&
&   'Likely your DDB is incorrect.',ch10,&
&   'Action: correct your DDB, or contact the ABINIT group.'
   MSG_ERROR(message)
 end if

end subroutine ddb_chkname
!!***

!----------------------------------------------------------------------

!!****f* m_ddb_hdr/compare_ddb_variables
!!
!! NAME
!! compare_ddb_variables
!!
!! FUNCTION
!! Compare the temporary DDB and input DDB preliminary information,
!! as well as psp information.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2019 ABINIT group (XG,MT,GA)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! NOTES
!! 1. All the variables have their usual meaning.
!! 2. Note that fullinit==0  means that the input DDB has been
!! initialized by a ground state input file. Some comparison are
!! then not required.
!! 3. All variables with 8 appended are from the new input DDB
!!
!! INPUTS
!!  acell, acell8 = lattice parameters
!!  amu, amu8 = atomic masses
!!  dimekb = dimension of KB projector set (only used for NCPP)
!!  ecut, ecut8 = cutoff energy
!!  ekb, ekb8 = KB energies for pseudopotentials
!!  fullinit, fullmrgddb_init = flags (see notes)
!!  iscf, iscf8 = SCF algorithm
!!  ixc, ixc8 = XC functional
!!  kpt, kpt8 = kpoint array
!!  kptnrm, kptnrm8 = normalization factor for kpt
!!  natom, natom8 = number of atoms
!!  nband, nband8 = number of bands at each kpt
!!  ngfft, ngfft8 = FFT grid sizes
!!  nkpt, nkpt8 = number of kpoints
!!  nsppol, nsppol8 = number of spin polarization (1 or 2)
!!  nsym, nsym8 = number of symmetry operations
!!  ntypat, ntypat8 = number of types of atoms
!!  occ, occ8 = occupation numbers
!!  occopt, occop8 = occupation style (metal, insulator, smearing...)
!!  pawecutdg,pawecutdg8= cutoff energy used for the fine "double grid" (PAW only)
!!  pawtab,pawtab8= PAW tabulated data (PAW dataset)
!!  rprim, rprim8 = primitive vectors of unit cell (cartesian coordinates)
!!  dfpt_sciss, dfpt_sciss8 = scissor correction (Ha)
!!  symrel, symrel8 = symmetry operations in reciprocal space
!!  tnons, tnons8 = translations associated to symrel
!!  tolwfr, tolwfr8 = tolerance on convergence of wavefunctions
!!  typat, typat8 = array of atom types
!!  usepaw = flag for utilization of PAW
!!  wtk, wtk8 = weights of kpoints
!!  xred, xred8 = reduced coordinates of atoms
!!  zion, zion8 = ionic charges of nuclei
!!
!! OUTPUT (corresponding values, checked and/or set)
!!  acell, amu, dimekb, ecut, ekb, fullinit, iscf, ixc, kpt, kptnrm,
!!  natom, nband, ngfft, nkpt, nsppol, nsym, ntypat, occ, occopt,
!!  rprim, dfpt_sciss, symrel, tnons, tolwfr, typat, usepaw, wtk, xred, zion
!!
!! PARENTS
!!      m_ddb_hdr
!!
!! CHILDREN
!!
!! SOURCE


subroutine compare_ddb_variables(&
& matom, matom8, mtypat, mtypat8, mkpt, mkpt8,&
& mband, mband8, msym, msym8,&
& acell,acell8,amu,amu8,dimekb,ecut,ecut8,ekb,ekb8,&
& fullinit,fullmrgddb_init,iscf,iscf8,ixc,ixc8,kpt,kpt8,&
& kptnrm,kptnrm8,&
& natom,natom8,nband,nband8,ngfft,ngfft8,nkpt,nkpt8,&
& nsppol,nsppol8,nsym,nsym8,ntypat,ntypat8,occ,occ8,&
& occopt,occop8,pawecutdg,pawecutdg8,pawtab,pawtab8,&
& rprim,rprim8,dfpt_sciss,dfpt_sciss8,symrel,symrel8,&
& tnons,tnons8,tolwfr,tolwfr8,typat,typat8,usepaw,wtk,wtk8,&
& xred,xred8,zion,zion8)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: matom, matom8, mtypat, mtypat8, mkpt, mkpt8
 integer,intent(in) :: mband, mband8, msym, msym8
 integer,intent(in) :: dimekb,fullmrgddb_init,iscf8,ixc8,natom8,nkpt8,nsppol8
 integer,intent(in) :: nsym8,ntypat8,occop8,usepaw
 integer,intent(inout) :: fullinit,iscf,ixc,natom,nkpt,nsppol,nsym,ntypat
 integer,intent(inout) :: occopt
 real(dp),intent(in) :: ecut8,kptnrm8,pawecutdg8,dfpt_sciss8,tolwfr8
 real(dp),intent(inout) :: ecut,kptnrm,pawecutdg,dfpt_sciss,tolwfr
!arrays
 integer,intent(in) :: nband8(mkpt8*nsppol8),ngfft8(18)
 integer,intent(inout) :: nband(mkpt*nsppol),ngfft(18)
 integer,intent(in) :: symrel8(3,3,msym8),typat8(matom8)
 integer,intent(inout) :: symrel(3,3,msym),typat(matom)
 real(dp),intent(in) :: acell8(3),amu8(mtypat8),ekb8(dimekb,mtypat8)
 real(dp),intent(inout) :: acell(3),amu(mtypat),ekb(dimekb,mtypat)
 real(dp),intent(in) :: kpt8(3,mkpt8),occ8(mband8*mkpt8*nsppol8)
 real(dp),intent(inout) :: kpt(3,mkpt),occ(mband*mkpt*nsppol)
 real(dp),intent(in) :: rprim8(3,3),tnons8(3,msym8)
 real(dp),intent(inout) :: rprim(3,3),tnons(3,msym)
 real(dp),intent(in) :: wtk8(mkpt8),xred8(3,matom8),zion8(mtypat8)
 real(dp),intent(inout) :: wtk(mkpt),xred(3,matom),zion(mtypat)
 type(pawtab_type),intent(in) :: pawtab8(ntypat8*usepaw)
 type(pawtab_type),intent(inout) :: pawtab(ntypat*usepaw)

!Local variables -------------------------
!scalars
 integer :: bantot,ii,ij,isym,itypat
 real(dp) :: ekbcm8,ekbcmp
 real(dp),parameter :: tol=2.0d-14
 character(len=500) :: msg

! *********************************************************************


!Compare all the preliminary information
!1. natom
 call chki8(natom,natom8,' natom')
!2. nkpt
!Compares the input and transfer values only if the input has not
!been initialized by a ground state input file
!There can also be the case of perturbation at Gamma, that
!only need half of the number of k points.
 if(fullinit/=0)then
   if(nkpt/=2*nkpt8 .and. 2*nkpt/=nkpt8)then

     ! GKA: We don't always need this variable to be consistent
     !      For example, we might have reduced the number of k-points
     !      with TRS only for certain q-points.
     !call chki8(nkpt,nkpt8,'  nkpt')
   else
     write(std_out,*)' compar8 : assume that one of the DDB to be',&
&     ' merged use Time-Reversal to'
     write(std_out,*)' decrease the number of k-points'
   end if
 else
!  Otherwise, takes the meaningful value
   nkpt=nkpt8
 end if
!3a. occopt
!Because the program will stop if the bloks
!do not compare well, take here the most favorable case.
 if(occop8==0)occopt=0
!3b. nband
!Compares the input and transfer values only if the input has not
!been initialized by a ground state input file
!There can also be the case of perturbation at Gamma, that
!only need half of the number of k points.
 if(fullinit==0 .or. nkpt8==2*nkpt)then
   bantot=0
   do ii=1,nkpt8
     nband(ii)=nband8(ii)
     bantot=bantot+nband(ii)
   end do
 else
   bantot=0
   do ii=1,nkpt
     if(nkpt==nkpt8)then
       call chki8(nband(ii),nband8(ii),' nband')
     end if
     bantot=bantot+nband(ii)
   end do
 end if
!9. nsppol
 call chki8(nsppol,nsppol8,'nsppol')
!4. nsym
 if(nsym/=1 .and. nsym8/=1)then
   call chki8(nsym,nsym8,'  nsym')
 end if
!5. ntypat
 call chki8(ntypat,ntypat8,'ntypat')
!6. acell
 do ii=1,3
   call chkr8(acell(ii),acell8(ii),' acell',tol)
 end do
!7. amu
 do ii=1,ntypat
   call chkr8(amu(ii),amu8(ii),'   amu',tol)
 end do
!9. date
!10. ecut
 call chkr8(ecut,ecut8,'  ecut',tol)
!10b. pawecutdg (PAW only)
 if (usepaw==1) then
   call chkr8(pawecutdg,pawecutdg8,'  ecut',tol)
 end if
!11. iscf
!Compares the input and transfer values only if the input has not
!been initialized by a ground state input file
 if(fullinit/=0)then
   call chki8(iscf,iscf8,'  iscf')
 else
!  Otherwise, takes the meaningful value
   iscf=iscf8
 end if
!12. ixc
 call chki8(ixc,ixc8,'   ixc')
!13. kpt and 14. kptnrm
!Compares the input and transfer values only if the input
!has not been initialized by a ground state input file
!and if the number of k points is identical
 if(nkpt8 == 2*nkpt .or. fullinit==0)then
!  Copy the largest number of k points in the right place
   do ij=1,nkpt8
     do ii=1,3
       kpt(ii,ij)=kpt8(ii,ij)
     end do
   end do
   kptnrm=kptnrm8
 else if (nkpt==nkpt8)then
   do ij=1,nkpt
     do ii=1,3
!      Compares the input and transfer values only if the input
!      has not been initialized by a ground state input file
       call chkr8(kpt(ii,ij)/kptnrm,kpt8(ii,ij)/kptnrm8,'   kpt',tol)
     end do
   end do
 end if
!16. ngfft
!MT dec 2013: deactivate the stop on ngfft to allow for
! (nfft-converged) DFPT calculations with GS WFK obtained with a different ngfft
 do ii=1,3
   if (ngfft(ii) == ngfft8(ii)) cycle
   write(msg,'(3a,i10,3a,i10,a)') &
&   'Comparing integers for variable ngfft.',ch10,&
&   'Value from input DDB is',ngfft(ii),' and',ch10,&
&   'from transfer DDB is',ngfft8(ii),'.'
   MSG_WARNING(msg)
 end do
!17. occ
!Compares the input and transfer values only if the input has not
!been inititialized by a ground state input file
 do ii=1,bantot
   if (fullinit==0 .or. nkpt8==2*nkpt) then
     occ(ii)=occ8(ii)
   else if(nkpt==nkpt8)then
     call chkr8(occ(ii),occ8(ii),'   occ',tol)
   end if
 end do
!18. rprim
 do ii=1,3
   do ij=1,3
     call chkr8(rprim(ii,ij),rprim8(ii,ij),' rprim',tol)
   end do
 end do
!19. dfpt_sciss
!Compares the input and transfer values only if the input has not
!been inititialized by a ground state input file
 if(fullinit/=0)then
   call chkr8(dfpt_sciss,dfpt_sciss8,' dfpt_sciss',tol)
 else
!  Otherwise, takes the meaningful value
   dfpt_sciss=dfpt_sciss8
 end if
!20. symrel
!If nsym == nsym8, compares the symmetry operations,
!otherwise, one of nsym or nsym8 is 1, and thus take the
!symrel corresponding to the largest set.
!nsym will be changed later
 if(nsym==nsym8)then
   do isym=1,nsym
     do ii=1,3
       do ij=1,3
         call chki8(symrel(ii,ij,isym),symrel8(ii,ij,isym),'symrel')
       end do
     end do
   end do
 else if(nsym8/=1)then
   symrel(:,:,1:nsym8)=symrel8(:,:,1:nsym8)
 end if
!21. tnons (see symrel)
 if(nsym==nsym8)then
   do isym=1,nsym
     do ii=1,3
       call chkr8(tnons(ii,isym),tnons8(ii,isym),' tnons',tol)
     end do
   end do
 else if(nsym8/=1)then
   tnons(:,1:nsym8)=tnons8(:,1:nsym8)
   nsym=nsym8
 end if
!22. tolwfr
!Take the less converged value...
 tolwfr=max(tolwfr,tolwfr8)
!23. typat
 do ii=1,ntypat
   call chki8(typat(ii),typat8(ii),' typat')
 end do
!24. wtk
!Compares the input and transfer values only if the input has not
!been initialized by a ground state input file and the
!number of k-points is identical.
 if(nkpt8==2*nkpt .or. fullinit==0)then
   do ii=1,nkpt8
     wtk(ii)=wtk8(ii)
   end do
 else if(nkpt==nkpt8)then
   do ii=1,nkpt
     call chkr8(wtk(ii),wtk8(ii),'   wtk',tol)
   end do
 end if
!25.xred
 do ij=1,natom
   do ii=1,3
     call chkr8(xred(ii,ij),xred8(ii,ij),'  xred',tol)
   end do
 end do
!26. zion
 do ii=1,ntypat
   call chkr8(zion(ii),zion8(ii),'  zion',tol)
 end do

!Finally, put the correct value of nkpt in the case
!of the use of the time-reversal symmetry
 if(2*nkpt==nkpt8)then
   nkpt=nkpt8
 end if

!Now compare the NC pseudopotential information
 if (usepaw==0) then
   if(dimekb/=0 .and. fullinit/=0 .and. fullmrgddb_init/=0 )then
     do ii=1,dimekb
       do itypat=1,ntypat
         ekbcmp=ekb(ii,itypat)
         ekbcm8=ekb8(ii,itypat)
         call chkr8(ekbcmp,ekbcm8,'   ekb',tol)
       end do
     end do
   else if(dimekb/=0 .and. fullmrgddb_init/=0)then
     do ii=1,dimekb
       do itypat=1,ntypat
         ekb(ii,itypat)=ekb8(ii,itypat)
       end do
     end do
   end if
 end if

!Now compare several PAW dataset information
 if (usepaw==1) then
   if (fullinit/=0 .and. fullmrgddb_init/=0) then
     do itypat=1,ntypat
       call chki8(pawtab(itypat)%basis_size,pawtab8(itypat)%basis_size,'bas_sz')
       call chki8(pawtab(itypat)%lmn_size,pawtab8(itypat)%lmn_size,'lmn_sz')
       call chki8(pawtab(itypat)%lmn2_size,pawtab8(itypat)%lmn2_size,'lmn2sz')
       call chkr8(pawtab(itypat)%rpaw,pawtab8(itypat)%rpaw,'  rpaw',tol3)
       call chkr8(pawtab(itypat)%rshp,pawtab8(itypat)%rshp,'rshape',tol3)
       call chki8(pawtab(itypat)%shape_type,pawtab8(itypat)%shape_type,'shp_tp')
       if (pawtab(itypat)%lmn2_size>0) then
         do ii=1,pawtab(itypat)%lmn2_size
           call chkr8(pawtab(itypat)%dij0(ii),pawtab8(itypat)%dij0(ii),'  dij0',tol)
         end do
       end if
     end do
   else if (fullmrgddb_init/=0) then
     do itypat=1,ntypat
       pawtab(itypat)%basis_size =pawtab8(itypat)%basis_size
       pawtab(itypat)%lmn_size   =pawtab8(itypat)%lmn_size
       pawtab(itypat)%rpaw       =pawtab8(itypat)%rpaw
       pawtab(itypat)%rshp       =pawtab8(itypat)%rshp
       pawtab(itypat)%shape_type =pawtab8(itypat)%shape_type
       if (pawtab8(itypat)%lmn2_size>0) then
         if (pawtab(itypat)%lmn2_size==0)  then
           ABI_ALLOCATE(pawtab(itypat)%dij0,(pawtab8(itypat)%lmn2_size))
         end if
         do ii=1,pawtab8(itypat)%lmn2_size
           pawtab(itypat)%dij0(ii)=pawtab8(itypat)%dij0(ii)
         end do
       end if
       pawtab(itypat)%lmn2_size  =pawtab8(itypat)%lmn2_size
     end do
   end if
 end if

end subroutine compare_ddb_variables
!!***

!----------------------------------------------------------------------

!!****f* m_ddb_hdr/chkr8
!!
!! NAME
!! chkr8
!!
!! FUNCTION
!! This small subroutine check the identity of reali and realt,
!! who are integers, and eventually send a message and stop
!! if they are found unequal by more than tol
!!
!! INPUTS
!! reali=first real number
!! intt=second  real number
!! character(len=6) name=name of the variable in the calling routine,
!!                       to be echoed
!! tol=tolerance
!!
!! OUTPUT
!!  (only checking)
!!
!! PARENTS
!!      m_ddb_hdr
!!
!! CHILDREN
!!
!! SOURCE

subroutine chkr8(reali,realt,name,tol)

!Arguments -------------------------------
!scalars
 real(dp),intent(in) :: reali,realt,tol
 character(len=6),intent(in) :: name

!Local variables-------------------------------
!scalars
 character(len=500) :: message

! *********************************************************************

 if(abs(reali-realt)>tol) then
   write(message, '(a,a,a,a,a,es16.6,a,a,a,es16.6,a,a,a)' )&
   'Comparing reals for variable',name,'.',ch10,&
   'Value from input DDB is',reali,' and',ch10,&
   'from transfer DDB is',realt,'.',ch10,&
   'Action: check your DDBs.'
   MSG_ERROR(message)
 end if

 end subroutine chkr8
!!***

!----------------------------------------------------------------------

!!****f* m_ddb_hdr/chki8
!!
!! NAME
!! chki8
!!
!! FUNCTION
!! This small subroutine check the identity of inti and intt,
!! who are integers, and eventually send a message and stop
!! if they are found unequal
!!
!! INPUTS
!! inti=first integer
!! intt=second integer
!! character(len=6) name=name of the variable in the calling routine,
!!                       to be echoed
!!
!! OUTPUT
!!  (only checking)
!!
!! PARENTS
!!      m_ddb_hdr
!!
!! CHILDREN
!!
!! SOURCE

subroutine chki8(inti,intt,name)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: inti,intt
 character(len=6),intent(in) :: name

!Local variables-------------------------------
!scalars
 character(len=500) :: message

! *********************************************************************

 if(inti/=intt) then
   write(message, '(a,a,a,a,a,i10,a,a,a,i10,a,a,a)' )&
   'Comparing integers for variable',name,'.',ch10,&
   'Value from input DDB is',inti,' and',ch10,&
   'from transfer DDB is',intt,'.',ch10,&
   'Action: check your DDBs.'
   MSG_ERROR(message)
 end if

 end subroutine chki8
!!***

!----------------------------------------------------------------------

!!****f* m_ddb_hdr/ddb_io_out
!!
!! NAME
!! ddb_io_out
!!
!! FUNCTION
!! Open Derivative DataBase, then
!! write Derivative DataBase preliminary information.
!! Note: only one processor writes the DDB.
!!
!! INPUTS
!! acell(3)=length scales of primitive translations (bohr)
!! amu(mtypat)=mass of the atoms (atomic mass unit)
!! dilatmx=the maximal dilatation factor
!! character(len=fnlen) dscrpt:string that describe the output database
!! ecut=kinetic energy planewave cutoff (hartree)
!! ecutsm=smearing energy for plane wave kinetic energy (Ha)
!! character(len=fnlen) filnam: name of output file
!! intxc=control xc quadrature
!! iscf=parameter controlling scf or non-scf choice
!! ixc=exchange-correlation choice parameter
!! kpt(3,mkpt)=k point set (reduced coordinates)
!! kptnrm=normalisation of k points
!! matom=maximum number of atoms
!! mband=maximum number of bands
!! mkpt=maximum number of special points
!! msym=maximum number of symetries
!! mtypat=maximum number of atom types
!! natom=number of atoms in the unit cell
!! nband(mkpt)=number of bands at each k point, for each polarization
!! ngfft(18)=contain all needed information about 3D FFT,
!!        see ~abinit/doc/variables/vargs.htm#ngfft
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
!! unddb=unit number for output
!! usepaw=flag for PAW
!! vrsddb=6 digit integer giving date, in form yymmdd for month=mm(1-12),
!!  day=dd(1-31), and year=yy(90-99 for 1990 to 1999,00-89 for 2000 to 2089),
!!  of current DDB version.
!! wtk(mkpt)=weight assigned to each k point
!! xred(3,matom)=reduced atomic coordinates
!! zion(mtypat)=valence charge of each type of atom
!! znucl(mtypat)=atomic number of atom type
!!
!! OUTPUT
!!  Only writing
!!
!! PARENTS
!!      m_ddb_hdr
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine ddb_io_out (dscrpt,filnam,matom,mband,&
&  mkpt,msym,mtypat,unddb,vrsddb,&
&  acell,amu,dilatmx,ecut,ecutsm,intxc,iscf,ixc,kpt,kptnrm,&
&  natom,nband,ngfft,nkpt,nspden,nspinor,nsppol,nsym,ntypat,occ,occopt,&
&  pawecutdg,rprim,dfpt_sciss,spinat,symafm,symrel,tnons,tolwfr,tphysel,tsmear,&
&  typat,usepaw,wtk,xred,zion,znucl)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: matom,mband,mkpt,msym,mtypat,unddb,vrsddb
 integer,intent(in) :: intxc,iscf,ixc,natom,nkpt,nspden,nspinor,nsppol,nsym
 integer,intent(in) :: ntypat,occopt,usepaw
 real(dp),intent(in) :: dilatmx,ecut,ecutsm,kptnrm,pawecutdg,dfpt_sciss,tolwfr,tphysel
 real(dp),intent(in) :: tsmear
 character(len=fnlen),intent(in) :: dscrpt,filnam
!arrays
 integer,intent(in) :: nband(mkpt*nsppol),ngfft(18),symafm(msym),symrel(3,3,msym)
 integer,intent(in) :: typat(matom)
 real(dp),intent(in) :: acell(3),amu(mtypat),kpt(3,mkpt),occ(mband*mkpt*nsppol)
 real(dp),intent(in) :: rprim(3,3),spinat(3,matom),tnons(3,msym),wtk(mkpt)
 real(dp),intent(in) :: xred(3,matom),zion(mtypat),znucl(mtypat)

!Local variables -------------------------
!Set routine version number here:
!scalars
 integer,parameter :: vrsio8=100401,vrsio8_old=010929,vrsio8_old_old=990527
 integer :: bantot,ii,ij,ikpt,iline,im,ierr
 character(len=500) :: message
!arrays
 character(len=9) :: name(9)

! *********************************************************************

 DBG_ENTER("COLL")


!Check ioddb8 version number (vrsio8) against mkddb version number
!(vrsddb)
 if (vrsio8/=vrsddb) then
   write(message, '(a,a,a,i10,a,a,i10,a)' )&
   ' ddb_io_out: WARNING -',ch10,&
   '  The input/output DDB version number=',vrsio8,ch10,&
   '  is not equal to the DDB version number=',vrsddb,'.'
   call wrtout(std_out,message,'COLL')
 end if

!Open the output derivative database.
!(version 2.1. : changed because of a bug in a Perl script
!should set up a name checking procedure, with change of name
!like for the output file)
 ierr = open_file(filnam,message,unit=unddb,status='unknown',form='formatted')
 if (ierr /= 0) then
   MSG_ERROR(message)
 end if

!Write the heading
 write(unddb, '(/,a,/,a,i10,/,/,a,a,/)' ) &
 ' **** DERIVATIVE DATABASE ****    ',&
 '+DDB, Version number',vrsddb,' ',trim(dscrpt)

!Write the descriptive data
!1. usepaw
 write(unddb, '(1x,a9,i10)' )'   usepaw',usepaw
!2. natom
 write(unddb, '(1x,a9,i10)' )'    natom',natom
!3. nkpt
 write(unddb, '(1x,a9,i10)' )'     nkpt',nkpt
!4. nsppol
 write(unddb, '(1x,a9,i10)' )'   nsppol',nsppol
!5. nsym
 write(unddb, '(1x,a9,i10)' )'     nsym',nsym
!6. ntypat
 write(unddb, '(1x,a9,i10)' )'   ntypat',ntypat
!7. occopt
 write(unddb, '(1x,a9,i10)' )'   occopt',occopt
!8. nband
 if(occopt==2)then
   im=12
   name(1)='    nband'
   do iline=1,(nkpt+11)/12
     if(iline==(nkpt+11)/12)im=nkpt-12*(iline-1)
     write(unddb, '(1x,a9,5x,12i5)' )name(1),(nband((iline-1)*12+ii),ii=1,im)
     name(1)='         '
   end do
   bantot=0
   do ikpt=1,nkpt
     bantot=bantot+nband(ikpt)
   end do
 else
   write(unddb, '(1x,a9,i10)' )'    nband',nband(1)
   bantot=nkpt*nband(1)
 end if

!9. acell
 write(unddb, '(1x,a9,3d22.14)' )'    acell',acell
!10. amu
 im=3
 name(1)='      amu'
 do iline=1,(ntypat+2)/3
   if(iline==(ntypat+2)/3)im=ntypat-3*(iline-1)
   write (unddb, '(1x,a9,3d22.14)' )name(1),(amu((iline-1)*3+ii),ii=1,im)
   name(1)='         '
 end do
!11. dilatmx
 write(unddb, '(1x,a9,d22.14)' )'  dilatmx',dilatmx
!12. ecut
 write(unddb, '(1x,a9,d22.14)' )'     ecut',ecut
!12b. pawecutdg (PAW)
 if (usepaw==1) then
   write(unddb, '(1x,a9,d22.14)' )'pawecutdg',pawecutdg
 end if
!13. ecutsm
 write(unddb, '(1x,a9,d22.14)' )'   ecutsm',ecutsm
!14. intxc
 write(unddb, '(1x,a9,i10)' )'    intxc',intxc
!15. iscf
 write(unddb, '(1x,a9,i10)' )'     iscf',iscf
!16. ixc
 write(unddb, '(1x,a9,i10)' )'      ixc',ixc
!17. kpt
 name(1)='      kpt'
 do iline=1,nkpt
   write (unddb, '(1x,a9,3d22.14)' )name(1),(kpt(ii,iline),ii=1,3)
   name(1)='      '
 end do
!18. kptnrm
 write(unddb, '(1x,a9,d22.14)' )'   kptnrm',kptnrm
!19. ngfft
 write(unddb, '(1x,a9,5x,3i5)' )'    ngfft',ngfft(1:3)
!20. nspden
 write(unddb, '(1x,a9,i10)' )'   nspden',nspden
!21. nspinor
 write(unddb, '(1x,a9,i10)' )'  nspinor',nspinor
!22. occ
 if(occopt==2)then
   im=3
   name(1)='      occ'
   do iline=1,(bantot+2)/3
     if(iline==(bantot+2)/3)im=bantot-3*(iline-1)
     write(unddb, '(1x,a9,3d22.14)' )name(1),(occ((iline-1)*3+ii),ii=1,im)
     name(1)='         '
   end do
 else
   im=3
   name(1)='      occ'
   do iline=1,(nband(1)+2)/3
     if(iline==(nband(1)+2)/3)im=nband(1)-3*(iline-1)
     write(unddb, '(1x,a9,3d22.14)' )name(1),(occ((iline-1)*3+ii),ii=1,im)
     name(1)='         '
   end do
 end if
!23. rprim
 name(1)='    rprim'
 do iline=1,3
   write(unddb, '(1x,a9,3d22.14)' )name(1),(rprim(ii,iline),ii=1,3)
   name(1)='      '
 end do
!24. dfpt_sciss
 write(unddb, '(1x,a11,d22.14)' )' dfpt_sciss',dfpt_sciss
!25. spinat
 name(1)='   spinat'
 do iline=1,natom
   write(unddb, '(1x,a9,3d22.14)' )name(1),(spinat(ii,iline),ii=1,3)
   name(1)='         '
 end do
!26. symafm
 im=12
 name(1)='   symafm'
 do iline=1,(nsym+11)/12
   if(iline==(nsym+11)/12)im=nsym-12*(iline-1)
   write(unddb, '(1x,a9,5x,12i5)' )name(1),(symafm((iline-1)*12+ii),ii=1,im)
   name(1)='         '
 end do
!27. symrel
 name(1)='   symrel'
 do iline=1,nsym
   write(unddb, '(1x,a9,5x,9i5)' )name(1),((symrel(ii,ij,iline),ii=1,3),ij=1,3)
   name(1)='         '
 end do
!28. tnons
 name(1)='    tnons'
 do iline=1,nsym
   write(unddb, '(1x,a9,3d22.14)' )name(1),(tnons(ii,iline),ii=1,3)
   name(1)='         '
 end do
!29. tolwfr
 write(unddb, '(1x,a9,d22.14)' )'   tolwfr',tolwfr
!30. tphysel
 write(unddb, '(1x,a9,d22.14)' )'  tphysel',tphysel
!31. tsmear
 write(unddb, '(1x,a9,d22.14)' )'   tsmear',tsmear
!32. typat
 im=12
 name(1)='    typat'
 do iline=1,(natom+11)/12
   if(iline==(natom+11)/12)im=natom-12*(iline-1)
   write(unddb, '(1x,a9,5x,12i5)' )name(1),(typat((iline-1)*12+ii),ii=1,im)
   name(1)='         '
 end do
!33. wtk
 name(1)='      wtk'
 im=3
 do iline=1,(nkpt+2)/3
   if(iline==(nkpt+2)/3)im=nkpt-3*(iline-1)
   write(unddb, '(1x,a9,3d22.14)' )name(1),(wtk((iline-1)*3+ii),ii=1,im)
   name(1)='         '
 end do
!34. xred
 name(1)='     xred'
 do iline=1,natom
   write(unddb, '(1x,a9,3d22.14)' )name(1),(xred(ii,iline),ii=1,3)
   name(1)='         '
 end do
!35. znucl
 name(1)='    znucl'
 im=3
 do iline=1,(ntypat+2)/3
   if(iline==(ntypat+2)/3)im=ntypat-3*(iline-1)
   write(unddb, '(1x,a9,3d22.14)' )name(1),(znucl((iline-1)*3+ii),ii=1,im)
   name(1)='         '
 end do
!36. zion
 name(1)='     zion'
 im=3
 do iline=1,(ntypat+2)/3
   if(iline==(ntypat+2)/3)im=ntypat-3*(iline-1)
   write(unddb, '(1x,a9,3d22.14)' )name(1),(zion((iline-1)*3+ii),ii=1,im)
   name(1)='         '
 end do

 DBG_EXIT("COLL")

end subroutine ddb_io_out
!!***

END MODULE m_ddb_hdr
!!***

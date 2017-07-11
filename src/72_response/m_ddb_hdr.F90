!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_ddb_hdr
!! NAME
!!  m_ddb
!!
!! FUNCTION
!!  This module contains the declaration of data types and methods
!!  to handle the header of the DDB files.
!!
!! COPYRIGHT
!! Copyright (C) 2011-2017 ABINIT group (GA)
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
 use m_xmpi

 use m_copy,    only : alloc_copy
 use m_pawtab,  only : pawtab_type, pawtab_nullify, pawtab_free !, pawtab_copy
 use m_psps,    only : psps_copy, psps_free
 use m_ddb,     only : psddb8, ioddb8_in, ddb_getdims

 implicit none

 private


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

subroutine ddb_hdr_init(ddb_hdr, dtset, psps, pawtab, ddb_version, dscrpt, &
                        nblok, xred, occ, ngfft)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ddb_hdr_init'
!End of the abilint section

 implicit none

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
   ddb_hdr%occ(:) = dtset%occ_orig(1:ddb_hdr%mband*ddb_hdr%mkpt*ddb_hdr%nsppol)
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

!----------------------------------------------------------------------

!!****f* m_ddb_hdr/ddb_hdr_malloc
!! NAME
!! ddb_hdr_malloc
!!
!! FUNCTION
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

subroutine ddb_hdr_malloc(ddb_hdr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ddb_hdr_malloc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(ddb_hdr_type),intent(inout) :: ddb_hdr

! ************************************************************************

 ! BEGIN DEBUG
 write(*,*) 'ddb_hdr_malloc: mkpt=', ddb_hdr%mkpt
 write(*,*) 'ddb_hdr_malloc: nsppol=', ddb_hdr%nsppol
 write(*,*) 'ddb_hdr_malloc: mtypat=', ddb_hdr%mtypat
 write(*,*) 'ddb_hdr_malloc: msym=', ddb_hdr%msym
 write(*,*) 'ddb_hdr_malloc: matom=', ddb_hdr%matom
 write(*,*) 'ddb_hdr_malloc: mband=', ddb_hdr%mband
 ! END DEBUG

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
!!
!! CHILDREN
!!
!! SOURCE

subroutine ddb_hdr_free(ddb_hdr)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ddb_hdr_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(ddb_hdr_type),intent(inout) :: ddb_hdr

! ************************************************************************

 ! integer
 if (allocated(ddb_hdr%nband)) then
   ABI_FREE(ddb_hdr%nband)
 end if
 if (allocated(ddb_hdr%symafm)) then
   ABI_FREE(ddb_hdr%symafm)
 end if
 if (allocated(ddb_hdr%symrel)) then
   ABI_FREE(ddb_hdr%symrel)
 end if
 if (allocated(ddb_hdr%typat)) then
   ABI_FREE(ddb_hdr%typat)
 end if

 ! real
 if (allocated(ddb_hdr%amu)) then
   ABI_FREE(ddb_hdr%amu)
 end if
 if (allocated(ddb_hdr%kpt)) then
   ABI_FREE(ddb_hdr%kpt)
 end if
 if (allocated(ddb_hdr%occ)) then
   ABI_FREE(ddb_hdr%occ)
 end if
 if (allocated(ddb_hdr%spinat)) then
   ABI_FREE(ddb_hdr%spinat)
 end if
 if (allocated(ddb_hdr%tnons)) then
   ABI_FREE(ddb_hdr%tnons)
 end if
 if (allocated(ddb_hdr%wtk)) then
   ABI_FREE(ddb_hdr%wtk)
 end if
 if (allocated(ddb_hdr%xred)) then
   ABI_FREE(ddb_hdr%xred)
 end if
 if (allocated(ddb_hdr%zion)) then
   ABI_FREE(ddb_hdr%zion)
 end if
 if (allocated(ddb_hdr%znucl)) then
   ABI_FREE(ddb_hdr%znucl)
 end if

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
!!
!! CHILDREN
!!
!! SOURCE

subroutine ddb_hdr_open_write(ddb_hdr, filnam, unddb, fullinit)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ddb_hdr_open_write'
 use interfaces_72_response
!End of the abilint section

 implicit none

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
!!
!! CHILDREN
!!
!! SOURCE

subroutine ddb_hdr_open_read(ddb_hdr, filnam, unddb, ddb_version, &
&        matom,mtypat,mband,mkpt,msym,dimekb,lmnmax,usepaw,dimonly)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ddb_hdr_open_read'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(ddb_hdr_type),intent(inout) :: ddb_hdr
 character(len=*),intent(in) :: filnam
 integer,intent(in) :: unddb
 integer,intent(in) :: ddb_version
 integer,intent(in),optional :: matom,mtypat,mband,mkpt,msym
 integer,intent(in),optional :: dimekb,lmnmax,usepaw
 integer,intent(in),optional :: dimonly

!Local variables -------------------------
 integer :: mblktyp,nblok
 integer :: matom_l,mtypat_l,mband_l,mkpt_l,msym_l
 integer :: dimekb_l,lmnmax_l,usepaw_l
 integer :: choice

! ************************************************************************

 ! Read the dimensions from the DDB
 call ddb_getdims(dimekb_l,filnam,lmnmax_l,mband_l,mblktyp, &
&       msym_l,matom_l,nblok,mkpt_l,mtypat_l,unddb,usepaw_l, &
&       ddb_version,xmpi_comm_self)

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

 ddb_hdr%mblktyp = mblktyp

 ddb_hdr%psps%dimekb = dimekb_l
 ddb_hdr%psps%ntypat = mtypat_l
 ddb_hdr%psps%lmnmax = lmnmax_l
 ddb_hdr%psps%usepaw = usepaw_l
 ddb_hdr%psps%useylm = usepaw_l  ! yep...

 if (present(dimonly) .and. (dimonly==1)) return

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
!!
!! CHILDREN
!!
!! SOURCE

subroutine ddb_hdr_compare(ddb_hdr1, ddb_hdr2)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ddb_hdr_compare'
 use interfaces_72_response
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(ddb_hdr_type),intent(inout) :: ddb_hdr1, ddb_hdr2

!Local variables -------------------------

! ************************************************************************

!    Should also compare indlmn and pspso ... but suppose that
!    the checking of ekb is enough for the psps.
!    Should also compare many other variables ... this is still
!    to be done ...
 call ddb_compare2 (&
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
END MODULE m_ddb_hdr


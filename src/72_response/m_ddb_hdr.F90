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

 use m_pawtab,  only : pawtab_type, pawtab_nullify, pawtab_copy, pawtab_free
 use m_psps,    only : psps_copy, psps_free
 use m_ddb,     only : psddb8

 implicit none

 private


 type,public :: ddb_hdr_type

   integer :: ddb_version
   ! Version of the DDB file

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
   ! nband(mkpt)

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
   ! occ(mband*mkpt)

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

subroutine ddb_hdr_init(ddb_hdr, dtset, psps, pawtab, ddb_version, &
&                       ngfft, occ, xred, dscrpt)


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
 integer,intent(in) :: ddb_version
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(in) :: xred(3,dtset%natom)
 character(len=fnlen),intent(in) :: dscrpt


! ************************************************************************

 ddb_hdr%ddb_version = ddb_version
 ddb_hdr%dscrpt = dscrpt
 ddb_hdr%ngfft = ngfft

 ! BEGIN DEBUG
 write(*,*) 'ddb_hdr_init: calling psps_copy'
 call flush()
 ! END DEBUG
 call psps_copy(psps, ddb_hdr%psps)

 ! BEGIN DEBUG
 write(*,*) 'ddb_hdr_init: will copy scalars'
 call flush()
 ! END DEBUG
 ! Copy scalars from dtset
 ddb_hdr%matom = dtset%natom
 ddb_hdr%natom = dtset%natom
 ddb_hdr%mband = dtset%mband
 ddb_hdr%mkpt = dtset%nkpt
 ddb_hdr%nkpt = dtset%nkpt
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

 call ddb_hdr_malloc(ddb_hdr)

 ! Copy arrays from dtset
 ddb_hdr%acell = dtset%acell_orig(1:3,1)
 ddb_hdr%rprim = dtset%rprim_orig(1:3,1:3,1)
 ddb_hdr%amu = dtset%amu_orig(:,1)
 ddb_hdr%nband = dtset%nband
 ddb_hdr%symafm = dtset%symafm
 ddb_hdr%symrel = dtset%symrel
 ddb_hdr%typat = dtset%typat
 ddb_hdr%kpt = dtset%kpt
 ddb_hdr%wtk = dtset%wtk
 ddb_hdr%spinat = dtset%spinat
 ddb_hdr%tnons = dtset%tnons
 ddb_hdr%zion = dtset%ziontypat
 ddb_hdr%znucl = dtset%znucl

 ddb_hdr%xred = xred
 ddb_hdr%occ = occ

 ! BEGIN DEBUG
 write(*,*) 'ddb_hdr_init: calling pawtab_copy'
 call flush()
 ! END DEBUG
 call pawtab_copy(pawtab, ddb_hdr%pawtab)
 ! BEGIN DEBUG
 write(*,*) 'ddb_hdr_init: done'
 call flush()
 ! END DEBUG

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

 ! integer
 ABI_MALLOC(ddb_hdr%nband,(ddb_hdr%mkpt))
 ABI_MALLOC(ddb_hdr%symafm,(ddb_hdr%msym))
 ABI_MALLOC(ddb_hdr%symrel,(3,3,ddb_hdr%msym))
 ABI_MALLOC(ddb_hdr%typat,(ddb_hdr%matom))

 ! real
 ABI_MALLOC(ddb_hdr%amu,(ddb_hdr%mtypat))
 ABI_MALLOC(ddb_hdr%kpt,(3, ddb_hdr%mkpt))
 ABI_MALLOC(ddb_hdr%occ,(ddb_hdr%mband*ddb_hdr%mkpt))
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
 ABI_FREE(ddb_hdr%nband)
 ABI_FREE(ddb_hdr%symafm)
 ABI_FREE(ddb_hdr%symrel)
 ABI_FREE(ddb_hdr%typat)

 ! real
 ABI_FREE(ddb_hdr%amu)
 ABI_FREE(ddb_hdr%kpt)
 ABI_FREE(ddb_hdr%occ)
 ABI_FREE(ddb_hdr%spinat)
 ABI_FREE(ddb_hdr%tnons)
 ABI_FREE(ddb_hdr%wtk)
 ABI_FREE(ddb_hdr%xred)
 ABI_FREE(ddb_hdr%zion)
 ABI_FREE(ddb_hdr%znucl)

 ! types
 call psps_free(ddb_hdr%psps)

 ! BEGIN DEBUG
 write(*,*) 'ddb_hdr_free: calling pawtab_free'
 call flush()
 ! END DEBUG
 if (allocated(ddb_hdr%pawtab)) then
   call pawtab_free(ddb_hdr%pawtab)
   ABI_DATATYPE_DEALLOCATE(ddb_hdr%pawtab)
 end if
 ! BEGIN DEBUG
 write(*,*) 'ddb_hdr_free: done'
 call flush()
 ! END DEBUG

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

subroutine ddb_hdr_open_write(ddb_hdr, filnam, unddb)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ddb_hdr_open_write'
 use interfaces_72_response
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(ddb_hdr_type),intent(inout) :: ddb_hdr
 character(len=fnlen),intent(in) :: filnam
 integer,intent(in) :: unddb

!Local variables -------------------------
 integer :: nblok, fullinit, choice

! ************************************************************************

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

 nblok=1 ; fullinit=1 ; choice=2
 call psddb8(choice,ddb_hdr%psps%dimekb,ddb_hdr%psps%ekb,fullinit,&
&  ddb_hdr%psps%indlmn,ddb_hdr%psps%lmnmax,nblok,ddb_hdr%ntypat,unddb,&
&  ddb_hdr%pawtab,ddb_hdr%psps%pspso,ddb_hdr%psps%usepaw,ddb_hdr%psps%useylm,&
&  ddb_hdr%ddb_version)


end subroutine ddb_hdr_open_write
!!***

!----------------------------------------------------------------------

END MODULE m_ddb_hdr


!!****m* ABINIT/m_anaddb_dataset
!! NAME
!!  m_anaddb_dataset
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2020 ABINIT group (XG,JCC,CL,MVeithen,XW,MJV)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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

module m_anaddb_dataset

 use defs_basis
 use m_abicore
 use m_errors

 use m_fstrings,  only : next_token, rmquotes, sjoin, inupper
 use m_symtk,     only : mati3det
 use m_parser,    only : intagm, chkvars_in_string, instrng
 use m_ddb,       only : DDB_QTOL

 implicit none

 private

 public :: anaddb_dataset_type
 public :: anaddb_dtset_free
 public :: anaddb_init
 public :: outvars_anaddb
 public :: invars9
!!***

!----------------------------------------------------------------------

!!****t* m_anaddb_dataset/anaddb_dataset_type
!! NAME
!! anaddb_dataset_type
!!
!! FUNCTION
!! The anaddb_dataset_type structured datatype
!! gather all the input variables for the anaddb code.
!!
!! SOURCE

 type anaddb_dataset_type

! Variables should be declared on separated lines in order to reduce the occurence of git conflicts.
! Since all these input variables are described in the anaddb_help.html
! file, they are not described in length here ...
! Integer
  integer :: alphon
  integer :: asr
  integer :: brav
  integer :: chneut
  integer :: dieflag
  integer :: dipdip
  integer :: dossum
  integer :: ep_scalprod
  integer :: eivec
  integer :: elaflag
  integer :: elphflag
  integer :: enunit
  integer :: gkk2write
  integer :: gkk_rptwrite
  integer :: gkqwrite
  integer :: gruns_nddbs
  integer :: ifcana
  integer :: ifcflag
  integer :: ifcout
  integer :: ifltransport
  integer :: instrflag
  integer :: natfix
  integer :: natifc
  integer :: natom
  integer :: natprj_bs
  integer :: nchan
  integer :: ndivsm=20
  integer :: nfreq
  integer :: ngrids
  integer :: nlflag
  integer :: nph1l
  integer :: nph2l
  integer :: nqpath
  integer :: nqshft
  integer :: nsphere
  integer :: nstrfix
  integer :: ntemper
  integer :: nwchan
  integer :: outboltztrap
  integer :: piezoflag
  integer :: polflag
  integer :: prtdos
  integer :: prt_ifc
  integer :: prtddb
  integer :: prtmbm
  integer :: prtfsurf
  integer :: prtnest
  integer :: prtphbands
  integer :: prtsrlr  ! print the short-range/long-range decomposition of phonon freq.
  integer :: prtvol = 0
  integer :: ramansr
  integer :: relaxat
  integer :: relaxstr
  integer :: rfmeth
  integer :: selectz
  integer :: symdynmat
  integer :: telphint
  integer :: thmflag
  integer :: qgrid_type
  integer :: ep_b_min
  integer :: ep_b_max
  integer :: ep_int_gkk
  integer :: ep_keepbands
  integer :: ep_nqpt
  integer :: ep_nspline
  integer :: ep_prt_yambo
  integer :: symgkq
  integer :: use_k_fine
  integer :: prtbltztrp

  integer :: ngqpt(9)             ! ngqpt(9) instead of ngqpt(3) is needed in wght9.f
  integer :: istrfix(6)
  integer :: ng2qpt(3)
  integer :: qrefine(3)
  integer :: kptrlatt(3,3)
  integer :: kptrlatt_fine(3,3)
  integer :: thermal_supercell(3,3)

! Real(dp)
  real(dp) :: a2fsmear
  real(dp) :: band_gap
  real(dp) :: dosdeltae
  real(dp) :: dossmear
  real(dp) :: dostol
  real(dp) :: elphsmear
  real(dp) :: elph_fermie
  real(dp) :: ep_extrael
  real(dp) :: freeze_displ
  real(dp) :: frmax
  real(dp) :: frmin
  real(dp) :: temperinc
  real(dp) :: tempermin
  real(dp) :: thmtol
  real(dp) :: mustar
  real(dp) :: rifcsph

  real(dp) :: q1shft(3,4)
  real(dp) :: q2shft(3)
  real(dp) :: targetpol(3)
  real(dp) :: vs_qrad_tolkms(2) = 0

! Integer arrays
  integer, allocatable :: atifc(:)
   ! atifc(natom) WARNING : there is a transformation of this input variable, in chkin9
   ! This should be changed ...

  integer, allocatable :: iatfix(:)
  ! iatfix(natom)

  integer, allocatable :: iatprj_bs(:)

! Real arrays
  real(dp), allocatable :: qnrml1(:)
  ! qnrml1(nph1l)

  real(dp), allocatable :: qnrml2(:)
  ! qnrml2(nph2l)

  real(dp), allocatable :: qpath(:,:)
  ! qpath(3,nqpath)

  real(dp), allocatable :: qph1l(:,:)
  ! qph1l(3,nph1l)

  real(dp), allocatable :: qph2l(:,:)
  ! qph2l(3,nph2l)

  real(dp), allocatable :: ep_qptlist(:,:)
  ! qph2l(3,ep_nqpt)

  character(len=fnlen), allocatable :: gruns_ddbs(:)
  ! gruns_ddbs(gruns_nddbs)

 end type anaddb_dataset_type
!!***

contains
!!***

!!****f* m_anaddb_dataset/anaddb_dtset_free
!!
!! NAME
!!   anaddb_dtset_free
!!
!! FUNCTION
!!   deallocate remaining arrays in the anaddb_dtset datastructure
!!
!! INPUTS
!!  anaddb_dtset = anaddb datastructure
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!      chkvars_in_string,inupper
!!
!! NOTES
!!
!! SOURCE

subroutine anaddb_dtset_free(anaddb_dtset)

!Arguments ------------------------------------
!scalars
 type(anaddb_dataset_type), intent(inout) :: anaddb_dtset

! *************************************************************************

 ABI_SFREE(anaddb_dtset%atifc)
 ABI_SFREE(anaddb_dtset%iatfix)
 ABI_SFREE(anaddb_dtset%iatprj_bs)
 ABI_SFREE(anaddb_dtset%qnrml1)
 ABI_SFREE(anaddb_dtset%qnrml2)
 ABI_SFREE(anaddb_dtset%qpath)
 ABI_SFREE(anaddb_dtset%qph1l)
 ABI_SFREE(anaddb_dtset%qph2l)
 ABI_SFREE(anaddb_dtset%ep_qptlist)
 ABI_SFREE(anaddb_dtset%gruns_ddbs)

end subroutine anaddb_dtset_free
!!***

!----------------------------------------------------------------------

!!****f* m_anaddb_dataset/invars9
!!
!! NAME
!! invars9
!!
!! FUNCTION
!! Open input file for the anaddb code, then reads or echoes the input information.
!!
!! INPUTS
!! lenstr=actual length of string
!! natom=number of atoms, needed for atifc
!! string*(*)=string of characters containing all input variables and data
!!
!! OUTPUT
!! anaddb_dtset= (derived datatype) contains all the input variables
!!
!! NOTES
!! Should be executed by one processor only.
!!
!! 27/01/2009: MJV: I have cleaned this routine extensively, putting all
!!  variables in alphabetical order, and in a second segment the dependent
!!  variables which need to be allocated depending on the dimensions read in.
!!  Could be divided into two routines as in abinit.
!!    FIXME: move checks to chkin9?
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!      chkvars_in_string,inupper
!!
!! SOURCE

subroutine invars9 (anaddb_dtset,lenstr,natom,string)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: lenstr,natom
 character(len=*),intent(in) :: string
 type(anaddb_dataset_type),intent(inout) :: anaddb_dtset

!Local variables -------------------------
!scalars
 integer,parameter :: vrsddb=100401 !Set routine version number here:
 integer,parameter :: jdtset=1
 integer :: ii,iph1,iph2,marr,tread,start
 integer :: idet
 character(len=500) :: message
 character(len=fnlen) :: path
!arrays
 integer,allocatable :: intarr(:)
 real(dp),allocatable :: dprarr(:)

!*********************************************************************
 marr=3
 ABI_ALLOCATE(intarr,(marr))
 ABI_ALLOCATE(dprarr,(marr))

!copy natom to anaddb_dtset
 anaddb_dtset%natom=natom

!=====================================================================
!start reading in dimensions and non-dependent variables
!=====================================================================

!A

!typical value for gaussian smearing of a2F function
 anaddb_dtset%a2fsmear = 0.00002_dp
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'a2fsmear',tread,'ENE')
 if(tread==1) anaddb_dtset%a2fsmear=dprarr(1)
 if (anaddb_dtset%a2fsmear < tol6) then
   write(message,'(a,f10.3,a,a,a,a,a)' )&
   'a2fsmear is ',anaddb_dtset%a2fsmear,', but only values > 1.e-6 ',ch10,&
   'are allowed',ch10,'Action: correct a2fsmear in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%alphon=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'alphon',tread,'INT')
 if(tread==1) anaddb_dtset%alphon=intarr(1)
!FIXME: need a test on input value

 anaddb_dtset%asr=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'asr',tread,'INT')
 if(tread==1) anaddb_dtset%asr=intarr(1)
 if(anaddb_dtset%asr<-2.or.anaddb_dtset%asr>5)then
   write(message, '(a,i0,5a)' )&
   'asr is ',anaddb_dtset%asr,', but the only allowed values',ch10,&
   'are 0, 1, 2, 3, 4, 5, -1 or -2 .',ch10,'Action: correct asr in your input file.'
!  Note : negative values are allowed when the acoustic sum rule
!  is to be applied after the analysis of IFCs
!  3,4 are for rotational invariance (under development)
!  5 is for hermitian imposition of the ASR
   MSG_ERROR(message)
 end if

!B

!Target band gap in eV
!The default value is just a very large number that will not be used in changing the band gap
 anaddb_dtset%band_gap = 999.0d0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'band_gap',tread,'DPR')
 if(tread==1) anaddb_dtset%band_gap=dprarr(1)

 anaddb_dtset%brav=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'brav',tread,'INT')
 if(tread==1) anaddb_dtset%brav=intarr(1)
 if(anaddb_dtset%brav<=-2.or.anaddb_dtset%brav>=5 .or. anaddb_dtset%brav==0)then
   write(message, '(a,i0,a5)' )&
   'brav is ',anaddb_dtset%brav,', but the only allowed values',ch10,&
   'are -1, 1,2,3 or 4 .',ch10,'Action: correct brav in your input file.'
   MSG_ERROR(message)
 end if

!C

 anaddb_dtset%chneut=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'chneut',tread,'INT')
 if(tread==1) anaddb_dtset%chneut=intarr(1)
 if(anaddb_dtset%chneut<0.or.anaddb_dtset%chneut>2)then
   write(message, '(a,i0,5a)' )&
   'chneut is ',anaddb_dtset%chneut,', but the only allowed values',ch10,&
   'are 0, 1 or 2.',ch10,'Action: correct chneut in your input file.'
   MSG_ERROR(message)
 end if

!D

 anaddb_dtset%dieflag=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dieflag',tread,'INT')
 if(tread==1) anaddb_dtset%dieflag=intarr(1)
 if(anaddb_dtset%dieflag<0.or.anaddb_dtset%dieflag>4)then
   write(message, '(a,i0,5a)' )&
   'dieflag is ',anaddb_dtset%dieflag,', but the only allowed values',ch10,&
   'are 0, 1, 2, 3 or 4.',ch10,'Action: correct dieflag in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%dipdip=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dipdip',tread,'INT')
 if(tread==1) anaddb_dtset%dipdip=intarr(1)
 if(anaddb_dtset%dipdip<0.or.anaddb_dtset%dipdip>1)then
   write(message, '(a,i0,5a)' )&
   'dipdip is ',anaddb_dtset%dipdip,', but the only allowed values',ch10,&
   'are 0 or 1 .',ch10,'Action: correct dipdip in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%ep_scalprod = 0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ep_scalprod',tread,'INT')
 if(tread==1) anaddb_dtset%ep_scalprod = intarr(1)
 if(anaddb_dtset%ep_scalprod < 0 .or. anaddb_dtset%ep_scalprod > 1) then
   write(message, '(a,i0,5a)' )&
   'ep_scalprod is ',anaddb_dtset%ep_scalprod,', but the only allowed values',ch10,&
   'are 0 or 1.',ch10,'Action: correct ep_scalprod in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%dosdeltae=one/Ha_cmm1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dosdeltae',tread,'DPR')
 if(tread==1) anaddb_dtset%dosdeltae=dprarr(1)

!FIXME : should probably be smaller
 anaddb_dtset%dossmear=5.0/Ha_cmm1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dossmear',tread,'DPR')
 if(tread==1) anaddb_dtset%dossmear=dprarr(1)
 if(anaddb_dtset%dossmear<=zero)then
   write(message, '(a,es14.4,3a)' )&
   'dossmear is ',anaddb_dtset%dossmear,', which is lower than 0 .',ch10,&
   'Action: correct dossmear in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%dostol=0.25_dp
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dostol',tread,'DPR')
 if(tread==1) anaddb_dtset%dostol=dprarr(1)
 if(anaddb_dtset%dostol<zero)then
   write(message, '(a,es14.4,3a)' )&
   'dostol is ',anaddb_dtset%dostol,', which is lower than 0 .',ch10,&
   'Action: correct dostol in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%dossum=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dossum',tread,'INT')
 if(tread==1) anaddb_dtset%dossum=intarr(1)
 if(anaddb_dtset%dossum < 0 .or. anaddb_dtset%dossum > one)then
   write(message, '(a,i0,5a)' )&
   'dossum is ',anaddb_dtset%dossum,', but the only allowed values',ch10,&
   'are 0, 1',ch10,'Action: correct dossum in your input file.'
   MSG_ERROR(message)
 end if

!E

 anaddb_dtset%eivec=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'eivec',tread,'INT')
 if(tread==1) anaddb_dtset%eivec=intarr(1)
 if(anaddb_dtset%eivec<0.or.anaddb_dtset%eivec>4)then
   write(message, '(a,i0,5a)' )&
   'eivec is ',anaddb_dtset%eivec,', but the only allowed values',ch10,&
   'are 0, 1, 2, 3 or 4.',ch10,'Action: correct eivec in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%elaflag=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'elaflag',tread,'INT')
 if(tread==1) anaddb_dtset%elaflag=intarr(1)
 if(anaddb_dtset%elaflag<0.or.anaddb_dtset%elaflag>5)then
   write(message,'(a,i0,5a)' )&
   'elaflag is ',anaddb_dtset%elaflag,', but the only allowed values',ch10,&
   'are 0,1,2,3,4 or 5 .',ch10,'Action: correct elaflag in your input file.'
   MSG_ERROR(message)
 end if

!By default use the real fermie (tests for abs(elph_fermie) < tol10 in the code)
 anaddb_dtset%elph_fermie = zero
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'elph_fermie',tread,'ENE')
 if(tread==1) anaddb_dtset%elph_fermie=dprarr(1)

!extra charge in unit cell (number of electrons) wrt neutral cell
!holes are negative values (reduce number of electrons)
 anaddb_dtset%ep_extrael = zero
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ep_extrael',tread,'DPR')
 if(tread==1) anaddb_dtset%ep_extrael=dprarr(1)

!number to control the spline interpolation in RTA
 anaddb_dtset%ep_nspline = 20
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ep_nspline',tread,'INT')
 if(tread==1) anaddb_dtset%ep_nspline=intarr(1)
 if(anaddb_dtset%ep_nspline < 0 .or. anaddb_dtset%ep_nspline > 1000) then
   write(message, '(a,i0,5a)' )&
   'ep_nspline is ',anaddb_dtset%ep_nspline,', but this should not be ',ch10,&
   'negative or too large .',ch10,'Action: correct ep_nspline in your input file.'
   MSG_ERROR(message)
 end if

!interpolate gkk or gamma. It should be better to interpolate gkk onto the
!k_phon, since the integration weights will be treated the same way
 anaddb_dtset%ep_int_gkk = 0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ep_int_gkk',tread,'INT')
 if(tread==1) anaddb_dtset%ep_int_gkk = intarr(1)
 if(anaddb_dtset%ep_int_gkk < 0 .or. anaddb_dtset%ep_int_gkk > 1) then
   write(message, '(a,i0,5a)' )&
   'ep_int_gkk is ',anaddb_dtset%ep_int_gkk,', but the only allowed values',ch10,&
   'are 0 or 1.',ch10,'Action: correct ep_int_gkk in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%elphflag=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'elphflag',tread,'INT')
 if(tread==1) anaddb_dtset%elphflag=intarr(1)
 if(anaddb_dtset%elphflag<0.or.anaddb_dtset%elphflag>1)then
   write(message,'(a,i0,5a)' )&
   'elphflag = ',anaddb_dtset%elphflag,', but the allowed values',ch10,&
   'are 0, or 1.',ch10,'Action: correct elphflag in your input file.'
   MSG_ERROR(message)
 end if

!typical value for gaussian smearing, but can vary sensibly with the metal
 anaddb_dtset%elphsmear = 0.01_dp
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'elphsmear',tread,'ENE')
 if(tread==1) anaddb_dtset%elphsmear=dprarr(1)
 if (anaddb_dtset%elphsmear < tol6) then
   write(message,'(a,f10.3,5a)' )&
   'elphsmear is ',anaddb_dtset%elphsmear,'. Only values > 1.e-6 ',ch10,&
   'are allowed',ch10,'Action: correct elphsmear in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%enunit=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'enunit',tread,'INT')
 if(tread==1) anaddb_dtset%enunit=intarr(1)
 if(anaddb_dtset%enunit<0.or.anaddb_dtset%enunit>2)then
   write(message, '(a,i0,5a)' )&
   'enunit is ',anaddb_dtset%enunit,', but the only allowed values',ch10,&
   'are 0, 1 or 2.',ch10,'Action: correct enunit in your input file.'
   MSG_ERROR(message)
 end if

!Default is 0 - not used unless telphint==2
 anaddb_dtset%ep_b_max = 0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ep_b_max',tread,'INT')
 if(tread==1) then
   anaddb_dtset%ep_b_max = intarr(1)
   if(anaddb_dtset%ep_b_max < 1) then
     write(message, '(a,i0,5a)' )&
     'ep_b_max is ',anaddb_dtset%ep_b_max,', but the only allowed values',ch10,&
     'are between 1 and nband.',ch10,'Action: correct ep_b_max in your input file.'
     MSG_ERROR(message)
   end if
 end if

!Default is 0 - not used unless telphint==2
 anaddb_dtset%ep_b_min = 0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ep_b_min',tread,'INT')
 if(tread==1) then
   anaddb_dtset%ep_b_min = intarr(1)
   if(anaddb_dtset%ep_b_min < 1) then
     write(message, '(a,i0,5a)' )&
     'ep_b_min is ',anaddb_dtset%ep_b_min,', but the only allowed values',ch10,&
     'are between 1 and nband.',ch10,'Action: correct ep_b_min in your input file.'
     MSG_ERROR(message)
   end if
 end if

 anaddb_dtset%ep_keepbands = 0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ep_keepbands',tread,'INT')
 if(tread==1) anaddb_dtset%ep_keepbands = intarr(1)
 if(anaddb_dtset%ep_keepbands < 0 .or. anaddb_dtset%ep_keepbands > 1) then
   write(message, '(a,i0,5a)' )&
   'ep_keepbands is ',anaddb_dtset%ep_keepbands,', but the only allowed values',ch10,&
   'are 0 or 1 .',ch10,'Action: correct ep_keepbands in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%ep_nqpt=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ep_nqpt',tread,'INT')
 if(tread==1) anaddb_dtset%ep_nqpt = intarr(1)
 if(anaddb_dtset%ep_nqpt < 0) then
   write(message, '(a,i0,5a)' )&
   'ep_nqpt is ',anaddb_dtset%ep_nqpt,', but the only allowed values',ch10,&
   'are > 0.',ch10,'Action: correct ep_nqpt in your input file.'
   MSG_ERROR(message)
 end if

 if (anaddb_dtset%ep_nqpt > 0) then
   ABI_ALLOCATE(anaddb_dtset%ep_qptlist,(3,anaddb_dtset%ep_nqpt))
   if(3*anaddb_dtset%ep_nqpt>marr)then
     marr=3*anaddb_dtset%ep_nqpt
     ABI_DEALLOCATE(intarr)
     ABI_DEALLOCATE(dprarr)
     ABI_ALLOCATE(intarr,(marr))
     ABI_ALLOCATE(dprarr,(marr))
   end if
   anaddb_dtset%ep_qptlist(:,:)=zero
   call intagm(dprarr,intarr,jdtset,marr,3*anaddb_dtset%ep_nqpt,string(1:lenstr),'ep_qptlist',tread,'DPR')
   if(tread==1) then
     anaddb_dtset%ep_qptlist(1:3,1:anaddb_dtset%ep_nqpt)=&
     reshape(dprarr(1:3*anaddb_dtset%ep_nqpt),(/3,anaddb_dtset%ep_nqpt/))
   else
     write(message,'(3a)')&
     'ep_nqpt is non zero but ep_qptlist is absent ',ch10,&
     'Action: specify ep_qptlist in your input file.'
     MSG_ERROR(message)
   end if
 end if


!F
 anaddb_dtset%freeze_displ = zero
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'freeze_displ',tread,'DPR')
 if(tread==1) anaddb_dtset%freeze_displ=dprarr(1)


 anaddb_dtset%frmax=ten
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'frmax',tread,'DPR')
 if(tread==1) anaddb_dtset%frmax=dprarr(1)
 if (anaddb_dtset%frmax < 0) then
   write(message,'(a,f10.3,5a)' )&
   'frmax is ',anaddb_dtset%frmax,'. Only values > 0 ',ch10,&
   'are allowed',ch10,'Action: correct frmax in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%frmin=zero
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'frmin',tread,'DPR')
 if(tread==1) anaddb_dtset%frmin=dprarr(1)
 if (anaddb_dtset%frmin < 0) then
   write(message,'(a,f10.3,5a)' )&
   'frmin is ',anaddb_dtset%frmin,'. Only values > 0 ',ch10,&
   'are allowed',ch10,'Action: correct frmin in your input file.'
   MSG_ERROR(message)
 end if

!G

 anaddb_dtset%gkk2write = 0
!call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gkk2write',&
!& tread,'INT')
!if(tread==1) anaddb_dtset%gkk2write = intarr(1)
!if(anaddb_dtset%gkk2write < 0 .or. anaddb_dtset%gkk2write > 1) then
!write(message, '(a,a,a,i8,a,a,a,a,a)' )&
!&  ' invars9 : ERROR -',ch10,&
!&  '  gkk2write is',anaddb_dtset%gkk2write,&
!&  ', but the only allowed values',ch10,&
!&  '  are 0 or 1 .',ch10,&
!&  '  Action: correct gkk2write in your input file.'
!MSG_ERROR(message)
!end if

 anaddb_dtset%gkk_rptwrite = 0
!call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gkk_rptwrite',&
!& tread,'INT')
!if(tread==1) anaddb_dtset%gkk_rptwrite = intarr(1)
!if(anaddb_dtset%gkk_rptwrite < 0 .or. anaddb_dtset%gkk_rptwrite > 1) then
!write(message, '(a,a,a,i8,a,a,a,a,a)' )&
!&  ' invars9 : ERROR -',ch10,&
!&  '  gkk_rptwrite is',anaddb_dtset%gkk_rptwrite,&
!&  ', but the only allowed values',ch10,&
!&  '  are 0 or 1 .',ch10,&
!&  '  Action: correct gkk_rptwrite in your input file.'
!MSG_ERROR(message)
!end if

 anaddb_dtset%gkqwrite = 0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gkqwrite',tread,'INT')
 if(tread==1) anaddb_dtset%gkqwrite = intarr(1)
 if(anaddb_dtset%gkqwrite < 0 .or. anaddb_dtset%gkqwrite > 1) then
   write(message, '(a,i0,5a)' )&
   'gkqwrite is ',anaddb_dtset%gkqwrite,', but the only allowed values',ch10,&
   'are 0 or 1.',ch10,'Action: correct gkqwrite in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%gruns_nddbs = 0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'gruns_nddbs',tread,'INT')
 if (tread==1) anaddb_dtset%gruns_nddbs = intarr(1)

 if (anaddb_dtset%gruns_nddbs /= 0) then
   ! Read list of DDB paths.
   ABI_MALLOC(anaddb_dtset%gruns_ddbs, (anaddb_dtset%gruns_nddbs))
   start = index(string, "GRUNS_DDBS") + len("GRUNS_DDBS") + 1
   do ii=1,anaddb_dtset%gruns_nddbs
     if (next_token(string, start, path) /= 0) then
       MSG_ERROR(sjoin("Cannot find DDB path in input string:", ch10, string(start:)))
     end if
     anaddb_dtset%gruns_ddbs(ii) = rmquotes(path)
   end do
 end if

!H

!I

 anaddb_dtset%ifcana=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ifcana',tread,'INT')
 if(tread==1) anaddb_dtset%ifcana=intarr(1)
 if(anaddb_dtset%ifcana<0.or.anaddb_dtset%ifcana>1)then
   write(message, '(a,i0,5a)' )&
   'ifcana is ',anaddb_dtset%ifcana,', but the only allowed values',ch10,&
   'are 0 or 1.',ch10,'Action: correct ifcana in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%ifcflag=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ifcflag',tread,'INT')
 if(tread==1) anaddb_dtset%ifcflag=intarr(1)
 if(anaddb_dtset%ifcflag<0.or.anaddb_dtset%ifcflag>1)then
   write(message, '(a,i0,5a)' )&
   'ifcflag is ',anaddb_dtset%ifcflag,', but the only allowed values',ch10,&
   'are 0 or 1.',ch10,'Action: correct ifcflag in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%ifcout=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ifcout',tread,'INT')
 if(tread==1) anaddb_dtset%ifcout=intarr(1)
 if(anaddb_dtset%ifcout<-1)then
   write(message, '(a,i0,3a)' )&
   'ifcout is ',anaddb_dtset%ifcout,', which is lower than -1.',ch10,&
   'Action: correct ifcout in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%ifltransport = 0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ifltransport',tread,'INT')
 if(tread==1) anaddb_dtset%ifltransport = intarr(1)
 if(anaddb_dtset%ifltransport < 0 .or. anaddb_dtset%ifltransport > 3) then
   write(message, '(a,i0,5a)' )&
   'ifltransport is ',anaddb_dtset%ifltransport,', but the only allowed values',ch10,&
   'are 0 or 1 or 2 or 3.',ch10,'Action: correct ifltransport in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%instrflag=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'instrflag',tread,'INT')
 if(tread==1) anaddb_dtset%instrflag=intarr(1)
 if(anaddb_dtset%instrflag<0.or.anaddb_dtset%instrflag>1)then
   write(message,'(a,i0,5a)' )&
   'instrflag is ',anaddb_dtset%instrflag,', but the only allowed values',ch10,&
   'are 0, 1.',ch10,'Action: correct instrflag in your input file.'
   MSG_ERROR(message)
 end if

!J

!K

 anaddb_dtset%kptrlatt = 0
!why this test on reading in kptrlatt?
 marr = 9
 ABI_DEALLOCATE(intarr)
 ABI_DEALLOCATE(dprarr)
 ABI_ALLOCATE(intarr,(marr))
 ABI_ALLOCATE(dprarr,(marr))
 call intagm(dprarr,intarr,jdtset,marr,9,string(1:lenstr),'kptrlatt',tread,'INT')
 if(tread==1)anaddb_dtset%kptrlatt(1:3,1:3)=reshape(intarr(1:9),(/3,3/))
!NOTE: no a priori way to test the validity of the integers in kptrlatt

 anaddb_dtset%kptrlatt_fine(:,:)=0
 marr = 9
 ABI_DEALLOCATE(intarr)
 ABI_DEALLOCATE(dprarr)
 ABI_ALLOCATE(intarr,(marr))
 ABI_ALLOCATE(dprarr,(marr))
 call intagm(dprarr,intarr,jdtset,marr,9,string(1:lenstr),'kptrlatt_fine',tread,'INT')
 if(tread==1)anaddb_dtset%kptrlatt_fine(1:3,1:3)=reshape(intarr(1:9),(/3,3/))


!L

!M

!typical value for mustar, but can vary sensibly with the metal
 anaddb_dtset%mustar = 0.1_dp
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'mustar',tread,'DPR')
 if(tread==1) anaddb_dtset%mustar=dprarr(1)
 if (anaddb_dtset%mustar < zero) then
   write(message,'(a,f10.3,5a)' )&
   'mustar is ',anaddb_dtset%mustar,', but only positive values',ch10,&
   'are allowed',ch10,'Action: correct mustar in your input file.'
   MSG_ERROR(message)
 end if

!N

 anaddb_dtset%natfix=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'natfix',tread,'INT')
 if(tread==1) anaddb_dtset%natfix=intarr(1)
 if(anaddb_dtset%natfix > natom)then
   write(message, '(a,i0,2a,i0,3a)' )&
   'natfix is ',anaddb_dtset%natfix,', which is larger than natom',' (=',natom,')',ch10,&
   'Action: correct natfix in your input file.'
   MSG_ERROR(message)
 end if

 if(anaddb_dtset%natfix < 0)then
   write(message, '(a,i0,3a)' )&
   'natfix is ',anaddb_dtset%natfix,', which is < 0',ch10,&
   'Action: correct natfix in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%natifc=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'natifc',tread,'INT')
 if(tread==1) anaddb_dtset%natifc=intarr(1)
 if(anaddb_dtset%natifc<0)then
   write(message, '(a,i0,3a)' )&
   'natifc is ',anaddb_dtset%natifc,', which is lower than 0 .',ch10,&
   'Action: correct natifc in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%natprj_bs=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'natprj_bs',tread,'INT')
 if(tread==1) anaddb_dtset%natprj_bs=intarr(1)
 if(anaddb_dtset%natprj_bs<0 .or. anaddb_dtset%natprj_bs > natom)then
   write(message, '(a,i0,a,i0,2a)' )&
   'natprj_bs is ',anaddb_dtset%natprj_bs,', but must be between 0 and natom = ',natom,ch10,&
   'Action: correct natprj_bs in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%nchan=800
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nchan',tread,'INT')
 if(tread==1) anaddb_dtset%nchan=intarr(1)
!FIXME: check this - it should probably be .ge. 1, not 0
 if(anaddb_dtset%nchan <0)then
   write(message, '(a,i0,3a)' )&
   'nchan is ',anaddb_dtset%nchan,', which is lower than 0 .',ch10,&
   'Action: correct nchan in your input file.'
   MSG_ERROR(message)
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ndivsm',tread,'INT')
 if(tread==1) anaddb_dtset%ndivsm=intarr(1)
 if(anaddb_dtset%ndivsm <=0)then
   write(message, '(a,i0,3a)' )&
   'ndivsm is ',anaddb_dtset%ndivsm,', which is <= 0 .',ch10,&
   'Action: correct ndivsm in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%nfreq=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nfreq',tread,'INT')
 if(tread==1) anaddb_dtset%nfreq=intarr(1)
 if(anaddb_dtset%nfreq<0)then
   write(message, '(a,i0,3a)' )&
   'nfreq is ',anaddb_dtset%nfreq,', which is lower than 0 .',ch10,&
   'Action: correct nfreq in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%ng2qpt(:)=0
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'ng2qpt',tread,'INT')
 if(tread==1) anaddb_dtset%ng2qpt(:)=intarr(1:3)
 do ii=1,3
   if(anaddb_dtset%ng2qpt(ii)<0)then
     write(message, '(a,i0,a,i0,3a,i0,a)' )&
     'ng2qpt(',ii,') is ',anaddb_dtset%ng2qpt(ii),', which is lower than 0 .',ch10,&
     'Action: correct ng2qpt(',ii,') in your input file.'
     MSG_ERROR(message)
   end if
 end do

 anaddb_dtset%ngqpt(:)=0
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'ngqpt',tread,'INT')
 if(tread==1) anaddb_dtset%ngqpt(1:3)=intarr(1:3)
 do ii=1,3
   if(anaddb_dtset%ngqpt(ii)<0)then
     write(message, '(a,i0,a,i0,3a,i0,a)' )&
     'ngqpt(',ii,') is ',anaddb_dtset%ngqpt(ii),', which is lower than 0 .',ch10,&
     'Action: correct ngqpt(',ii,') in your input file.'
     MSG_ERROR(message)
   end if
 end do

 anaddb_dtset%ngrids=4
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ngrids',tread,'INT')
 if(tread==1) anaddb_dtset%ngrids=intarr(1)
 if(anaddb_dtset%ngrids<0)then
   write(message, '(a,i0,3a)' )&
   'ngrids is ',anaddb_dtset%ngrids,', which is lower than 0 .',ch10,&
   'Action: correct ngrids in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%nlflag=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nlflag',tread,'INT')
 if(tread==1) anaddb_dtset%nlflag=intarr(1)
 if(anaddb_dtset%nlflag<0.or.anaddb_dtset%nlflag>3)then
   write(message, '(a,i0,5a)' )&
   'nlflag is ',anaddb_dtset%nlflag,', but the only allowed values',ch10,&
   'are 0, 1, 2 or 3.',ch10,'Action: correct nlflag in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%nph1l=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nph1l',tread,'INT')
 if(tread==1) anaddb_dtset%nph1l=intarr(1)
 if(anaddb_dtset%nph1l<0)then
   write(message, '(a,i0,3a)' )&
   'nph1l is ',anaddb_dtset%nph1l,', which is lower than 0 .',ch10,&
   'Action: correct nph1l in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%nph2l=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nph2l',tread,'INT')
 if(tread==1) anaddb_dtset%nph2l=intarr(1)
 if(anaddb_dtset%nph2l<0)then
   write(message, '(a,i0,3a)' )&
   'nph2l is ',anaddb_dtset%nph2l,', which is lower than 0 .',ch10,&
   'Action: correct nph2l in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%nqpath=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nqpath',tread,'INT')
 if(tread==1) anaddb_dtset%nqpath=intarr(1)
 if(anaddb_dtset%nqpath<0)then
   write(message,'(a,i0,3a)' )&
   'nqpath is ',anaddb_dtset%nqpath,', but must be positive',ch10,&
   'Action: correct elphflag in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%nqshft=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nqshft',tread,'INT')
 if(tread==1) anaddb_dtset%nqshft=intarr(1)
 if(anaddb_dtset%nqshft<0 .or. anaddb_dtset%nqshft==3 .or. anaddb_dtset%nqshft>=5 )then
   write(message, '(a,i0,5a)' )&
   'nqshft is ',anaddb_dtset%nqshft,', but the only allowed values',ch10,&
   'are 1, 2 or 4 .',ch10,'Action: correct nqshft in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%nsphere=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nsphere',tread,'INT')
 if(tread==1) anaddb_dtset%nsphere=intarr(1)
 if(anaddb_dtset%nsphere < -1)then
   write(message, '(a,i0,3a)' )&
   'nsphere is ',anaddb_dtset%nsphere,', while it must be >= 0 or equal to -1',ch10,&
   'Action: correct nsphere in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%nstrfix=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nstrfix',tread,'INT')
 if(tread==1) anaddb_dtset%nstrfix=intarr(1)
 if(anaddb_dtset%nstrfix > 6)then
   write(message, '(a,i0,3a)' )&
   'nstrfix is ',anaddb_dtset%nstrfix,', which is larger than 6',ch10,&
   'Action: correct nstrfix in your input file.'
   MSG_ERROR(message)
 end if

 if(anaddb_dtset%nstrfix < 0)then
   write(message, '(a,i0,3a)' )&
   'nstrfix is ',anaddb_dtset%nstrfix,', which is < 0',ch10,&
   'Action: correct nstrfix in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%ntemper=10
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ntemper',tread,'INT')
 if(tread==1) anaddb_dtset%ntemper=intarr(1)
 if(anaddb_dtset%ntemper <0)then
   write(message, '(a,i0,3a)' )&
   'ntemper is ',anaddb_dtset%ntemper,', which is lower than 0',ch10,&
   'Action: correct ntemper in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%nwchan=10
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nwchan',tread,'INT')
 if(tread==1) anaddb_dtset%nwchan=intarr(1)
!FIXME: check this - it should probably be .ge. 1, not 0
 if(anaddb_dtset%nwchan<0)then
   write(message, '(a,i0,3a)' )&
   'nwchan is ',anaddb_dtset%nwchan,', which is lower than 0 .',ch10,&
   'Action: correct nwchan in your input file.'
   MSG_ERROR(message)
 end if

!O
 anaddb_dtset%outboltztrap = 0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'outboltztrap',tread,'INT')
 if(tread==1) anaddb_dtset%outboltztrap=intarr(1)
 if(anaddb_dtset%outboltztrap<0.or.anaddb_dtset%outboltztrap>1)then
   write(message,'(a,i0,5a)' )&
   'outboltztrap is ',anaddb_dtset%outboltztrap,', but the only allowed values',ch10,&
   'are 0 or 1.',ch10,'Action: correct outboltztrap in your input file.'
   MSG_ERROR(message)
 end if


!P
 anaddb_dtset%piezoflag=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'piezoflag',tread,'INT')
 if(tread==1) anaddb_dtset%piezoflag=intarr(1)
 if(anaddb_dtset%piezoflag<0.or.anaddb_dtset%piezoflag>7)then
   write(message,'(3a,i0,5a)' )&
   ' piezoflag is ',anaddb_dtset%piezoflag,', but the only allowed values',ch10,&
   'are 0, 1,2,3,4,5,6,7  .',ch10,'Action: correct piezoflag in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%polflag=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'polflag',tread,'INT')
 if(tread==1) anaddb_dtset%polflag=intarr(1)
 if(anaddb_dtset%polflag<0.or.anaddb_dtset%polflag>1)then
   write(message, '(a,i0,5a)' )&
   'polflag is ',anaddb_dtset%polflag,', but the only allowed values',ch10,&
   'are 0 or 1.',ch10,'Action: correct polflag in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%prtbltztrp=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtbltztrp',tread,'INT')
 if(tread==1) anaddb_dtset%prtbltztrp=intarr(1)
 if(anaddb_dtset%prtbltztrp<0.or.anaddb_dtset%prtbltztrp>1)then
   write(message, '(a,i0,5a)' )&
   'prtbltztrp is ',anaddb_dtset%prtbltztrp,', but the only allowed values',ch10,&
   'are 0 or 1.',ch10,'Action: correct prtbltztrp in your input file.'
   MSG_ERROR(message)
 end if

 ! Default is no output for PHDOS
 anaddb_dtset%prtdos=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtdos',tread,'INT')
 if(tread==1) anaddb_dtset%prtdos = intarr(1)
 if(anaddb_dtset%prtdos < 0 .or. anaddb_dtset%prtdos > 2) then
   write(message, '(a,i0,5a)' )&
   'prtdos is ',anaddb_dtset%prtdos,', but the only allowed values',ch10,&
   'are 0 (no output) or 1 (gaussians) or 2 (tetrahedra) ',ch10,&
   'Action: correct prtdos in your input file.'
   MSG_ERROR(message)
 end if

!Default is no output for the Fermi Surface
 anaddb_dtset%prtfsurf = 0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtfsurf',tread,'INT')
 if(tread==1) anaddb_dtset%prtfsurf = intarr(1)
 if(anaddb_dtset%prtfsurf < 0 .or. anaddb_dtset%prtfsurf > 2) then
   write(message, '(a,i0,5a)' )&
   'prtfsurf is ',anaddb_dtset%prtfsurf,'. The only allowed values',ch10,&
   'are 0 (no output) or 1 (Xcrysden bxsf format)',ch10,  &
   'Action: correct prtfsurf in your input file.'
   MSG_ERROR(message)
 end if

!Default is no output of the real space IFC to file
 anaddb_dtset%prt_ifc = 0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prt_ifc',tread,'INT')
 if(tread==1) anaddb_dtset%prt_ifc = intarr(1)
 if(anaddb_dtset%prt_ifc < 0 .or. anaddb_dtset%prt_ifc > 1) then
   write(message, '(a,i0,5a)' )&
   'prtf_ifc is ',anaddb_dtset%prt_ifc,'. The only allowed values',ch10,&
   'are 0 (no output) or 1 (AI2PS format)',ch10,  &
   'Action: correct prt_ifc in your input file.'
   MSG_ERROR(message)
 end if
! check that ifcout is set
 if (anaddb_dtset%prt_ifc /= 0 .and. anaddb_dtset%ifcout == 0) then
   anaddb_dtset%ifcout = -1 ! this forces output of all IFC
 end if

!Default is no output of the DDB to file
 anaddb_dtset%prtddb = 0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtddb',tread,'INT')
 if(tread==1) anaddb_dtset%prtddb = intarr(1)
 if(anaddb_dtset%prtddb < 0 .or. anaddb_dtset%prtddb > 1) then
   write(message, '(a,i0,5a)' )&
   'prtf_ddb is ',anaddb_dtset%prtddb,'. The only allowed values',ch10,&
   'are 0 (no output) or 1 (print DDB and DDB.nc files)',ch10,  &
   'Action: correct prtddb in your input file.'
   MSG_ERROR(message)
 end if
! check that ifcflag is set
 if (anaddb_dtset%prtddb /= 0 .and. anaddb_dtset%ifcflag == 0) then
   anaddb_dtset%ifcflag = 1 ! this forces the use of IFC
 end if

 anaddb_dtset%prtmbm=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtmbm',tread,'INT')
 if(tread==1) anaddb_dtset%prtmbm=intarr(1)
!FIXME: should check whether value of prtmbm is valid

!Default is no output of the nesting factor
 anaddb_dtset%prtnest = 0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtnest',tread,'INT')
 if(tread==1) anaddb_dtset%prtnest = intarr(1)
 if(anaddb_dtset%prtnest < 0 .or. anaddb_dtset%prtnest > 2) then
   write(message, '(a,i0,5a)' )&
   'prtnest is ',anaddb_dtset%prtnest,' The only allowed values',ch10,&
   'are 0 (no nesting), 1 (XY format) or 2 (XY + Xcrysden format)',ch10,&
   'Action: correct prtnest in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%prtphbands=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtphbands',tread,'INT')
 if (tread==1) anaddb_dtset%prtphbands=intarr(1)
 if (all(anaddb_dtset%prtphbands /= [0,1,2])) then
   write(message, '(a,i0,a)' )&
    'prtphbands is ',anaddb_dtset%prtphbands,', but the only allowed values are [0, 1, 2].'
   MSG_ERROR(message)
 end if

 anaddb_dtset%prtsrlr=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtsrlr',tread,'INT')
 if(tread==1) anaddb_dtset%prtsrlr=intarr(1)
 if(anaddb_dtset%prtsrlr<0.or.anaddb_dtset%prtsrlr>1)then
   write(message, '(a,i0,5a)' )&
   'prtsrlr is ',anaddb_dtset%prtsrlr,', but the only allowed values',ch10,&
   'are 0 or 1.',ch10,'Action: correct prtsrlr in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%prtvol = 0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtvol',tread,'INT')
 if(tread==1) anaddb_dtset%prtvol = intarr(1)

!Q

 anaddb_dtset%q2shft(:)=zero
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'q2shft',tread,'DPR')
 if(tread==1) anaddb_dtset%q2shft(:)=dprarr(1:3)
!FIXME: need a test on valid entries for q2shft

 anaddb_dtset%qgrid_type=1 ! default is uniform nqpt(:) grid
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'qgrid_type',tread,'INT')
 if(tread==1) anaddb_dtset%qgrid_type = intarr(1)
 if(anaddb_dtset%qgrid_type < 1 .or. anaddb_dtset%qgrid_type > 2) then
   write(message, '(a,i0,5a)' )&
   'qgrid_type is ',anaddb_dtset%qgrid_type,' The only allowed values',ch10,&
   'are 1 (uniform grid from nqpt) or 2 (listed in ep_nqpt, ep_qptlist)',ch10,&
   'Action: correct qgrid_type in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%qrefine=1 ! default is no refinement
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'qrefine',tread,'INT')
 if(tread==1) anaddb_dtset%qrefine = intarr(1:3)
 do ii=1,3
   if(anaddb_dtset%qrefine(ii) < 1) then
     write(message, '(a,3i0,a,a,a,a,a)' )&
     'qrefine is',anaddb_dtset%qrefine,' The only allowed values',ch10,&
     'are integers >= 1 giving the refinement of the ngqpt grid',ch10,&
     'Action: correct qrefine in your input file.'
     MSG_ERROR(message)
   end if
 end do

!R

 anaddb_dtset%ramansr=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ramansr',tread,'INT')
 if(tread==1) anaddb_dtset%ramansr=intarr(1)
 if(anaddb_dtset%ramansr<0.or.anaddb_dtset%ramansr>2)then
   write(message, '(a,i0,5a)' )&
   'ramansr is ',anaddb_dtset%ramansr,', but the only allowed values',ch10,&
   'are 0, 1 or 2.',ch10,'Action: correct ramansr in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%relaxat=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'relaxat',tread,'INT')
 if(tread==1) anaddb_dtset%relaxat=intarr(1)
 if(anaddb_dtset%relaxat < 0.or.anaddb_dtset%relaxat > 1)then
   write(message, '(a,i0,5a)' )&
   'relaxat is ',anaddb_dtset%relaxat,', but the only allowed values',ch10,&
   'are 0 or 1.',ch10,'Action: correct relaxat in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%relaxstr=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'relaxstr',tread,'INT')
 if(tread==1) anaddb_dtset%relaxstr=intarr(1)
 if(anaddb_dtset%relaxstr<0.or.anaddb_dtset%relaxstr>1)then
   write(message, '(a,i0,5a)' )&
   'relaxstr is ',anaddb_dtset%relaxstr,'but the only allowed values',ch10,&
   'are 0 or 1.',ch10,'Action: correct relaxstr in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%rfmeth=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rfmeth',tread,'INT')
 if(tread==1) anaddb_dtset%rfmeth=intarr(1)
 if(anaddb_dtset%rfmeth<1.or.anaddb_dtset%rfmeth>2)then
   write(message, '(a,i0,5a)' )&
   'rfmeth is ',anaddb_dtset%rfmeth,', but the only allowed values',ch10,&
   'are 1 or 2.',ch10,'Action: correct rfmeth in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%rifcsph=zero
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rifcsph',tread,'DPR')
 if(tread==1) anaddb_dtset%rifcsph=dprarr(1)
! if(anaddb_dtset%rifcsph<-tol12)then
!   write(message, '(a,f10.3,3a)' )&
!&   'rifcsph is ',anaddb_dtset%rifcsph,', which is lower than zero.',ch10,&
!&   'Action: correct rifcsph in your input file.'
!   MSG_ERROR(message)
! end if

!S

 anaddb_dtset%selectz=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'selectz',tread,'INT')
 if(tread==1) anaddb_dtset%selectz=intarr(1)
 if(anaddb_dtset%selectz<0.or.anaddb_dtset%selectz>2)then
   write(message, '(a,i0,5a)' )&
   'selectz is ',anaddb_dtset%selectz,', but the only allowed values',ch10,&
   'are 0, 1 or 2 .',ch10,'Action: correct selectz in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%symdynmat=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'symdynmat',tread,'INT')
 if(tread==1) anaddb_dtset%symdynmat=intarr(1)
 if(anaddb_dtset%symdynmat/=0.and.anaddb_dtset%symdynmat/=1)then
   write(message, '(a,i0,5a)' )&
   'symdynmat is ',anaddb_dtset%symdynmat,'. The only allowed values',ch10,&
   'are 0, or 1.',ch10,'Action: correct symdynmat in your input file.'
   MSG_ERROR(message)
 end if

!T

 anaddb_dtset%targetpol(:) = 0._dp
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'targetpol',tread,'DPR')
 if(tread==1) anaddb_dtset%targetpol(1:3) = dprarr(1:3)

!Default is use gaussian integration
 anaddb_dtset%telphint = 1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'telphint',tread,'INT')
 if(tread==1) anaddb_dtset%telphint = intarr(1)
 if(anaddb_dtset%telphint < 0 .or. anaddb_dtset%telphint > 3) then
   write(message, '(a,i0,6a)' )&
   'telphint is ',anaddb_dtset%telphint,'. The only allowed values',ch10,&
   'are 0 (tetrahedron) or 1 (gaussian) or ','2 (set of bands occupied ep_b_min,ep_b_max) or 3 (Fermi Dirac).',ch10,&
   'Action: correct telphint in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%temperinc=100.0_dp
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'temperinc',tread,'DPR')
 if(tread==1) anaddb_dtset%temperinc=dprarr(1)
 if(anaddb_dtset%temperinc < zero)then
   write(message, '(a,f10.3,3a)' )&
   'temperinc is ',anaddb_dtset%temperinc,', which is lower than 0 .',ch10,&
   'Action: correct temperinc in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%tempermin=100.0_dp
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'tempermin',tread,'DPR')
 if(tread==1) anaddb_dtset%tempermin=dprarr(1)
 if(anaddb_dtset%tempermin<-tol12)then
   write(message, '(a,f10.3,3a)' )&
   'tempermin is ',anaddb_dtset%tempermin,', which is lower than 0 .',ch10,&
   'Action: correct tempermin in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%thermal_supercell(:,:)=0
 marr = 9
 ABI_DEALLOCATE(intarr)
 ABI_DEALLOCATE(dprarr)
 ABI_ALLOCATE(intarr,(marr))
 ABI_ALLOCATE(dprarr,(marr))
 call intagm(dprarr,intarr,jdtset,marr,9,string(1:lenstr),'thermal_supercell',tread,'INT')
 if(tread==1) anaddb_dtset%thermal_supercell(1:3,1:3)=reshape(intarr(1:9),(/3,3/))
 call mati3det(anaddb_dtset%thermal_supercell, idet)
 if(sum(abs(anaddb_dtset%thermal_supercell))>0 .and. idet == 0) then
   write(message, '(a,9I6,5a)' )&
   'thermal_supercell is ',anaddb_dtset%thermal_supercell,', but the matrix must be non singular',ch10,&
   'with a non zero determinant.',ch10,'Action: correct thermal_supercell in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%thmflag=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'thmflag',tread,'INT')
 if(tread==1) anaddb_dtset%thmflag=intarr(1)
 if(anaddb_dtset%thmflag<0.or.anaddb_dtset%thmflag>8)then
   write(message, '(a,i0,5a)' )&
   'thmflag is ',anaddb_dtset%thmflag,', but the only allowed values',ch10,&
   'are between 0 to 8 (included).',ch10,'Action: correct thmflag in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%thmtol=0.25_dp
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'thmtol',tread,'DPR')
 if(tread==1) anaddb_dtset%thmtol=dprarr(1)
 if(anaddb_dtset%thmtol<zero)then
   write(message, '(a,es14.4,3a)' )&
   'thmtol is ',anaddb_dtset%thmtol,', which is lower than 0 .',ch10,&
   'Action: correct thmtol in your input file.'
   MSG_ERROR(message)
 end if

 anaddb_dtset%ep_prt_yambo = 0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ep_prt_yambo',tread,'INT')
 if(tread==1) anaddb_dtset%ep_prt_yambo = intarr(1)
 if(anaddb_dtset%ep_prt_yambo< 0 .or. anaddb_dtset%ep_prt_yambo> 1) then
   write(message, '(a,i0,5a)' )&
   'ep_prt_yambo is ',anaddb_dtset%ep_prt_yambo,', but the only allowed values',ch10,&
   'are 0 or 1.',ch10,'Action: correct ep_prt_yambo in your input file.'
   MSG_ERROR(message)
 end if

!default means _do_ symmetrize the ep coupling matrices over qpoints
 anaddb_dtset%symgkq = 1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'symgkq',tread,'INT')
 if(tread==1) anaddb_dtset%symgkq = intarr(1)
 if(anaddb_dtset%symgkq< 0 .or. anaddb_dtset%symgkq> 1) then
   write(message, '(a,i0,5a)' )&
   'symgkq is ',anaddb_dtset%symgkq,', but the only allowed values',ch10,&
   'are 0 or 1.',ch10,'Action: correct symgkq in your input file.'
   MSG_ERROR(message)
 else if (anaddb_dtset%symgkq == 0) then
   MSG_WARNING('You have turned off el-ph matrix symmetrization over q. Use at own risk')
 end if

!U

 anaddb_dtset%use_k_fine = 0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'use_k_fine',tread,'INT')
 if(tread==1) anaddb_dtset%use_k_fine= intarr(1)
 if(anaddb_dtset%use_k_fine /= 1 .and. anaddb_dtset%use_k_fine /= 0) then
   write(message, '(a,i0,5a)' )&
   'use_k_fine is ',anaddb_dtset%use_k_fine,', but the only allowed values',ch10,&
   'are 1 or 0.',ch10,'Action: correct use_k_fine in your input file.'
   MSG_ERROR(message)
 end if

 if(anaddb_dtset%use_k_fine == 1) then
   if (sum(anaddb_dtset%kptrlatt) == 0 .or. sum(anaddb_dtset%kptrlatt_fine) == 0 ) then
     MSG_ERROR('If a finer k-grid is used, you must specify both kptrlatt and kptrlatt_fine')
   end if
 end if


!V
 anaddb_dtset%vs_qrad_tolkms(:) = zero
 call intagm(dprarr,intarr,jdtset,marr,2,string(1:lenstr),'vs_qrad_tolkms',tread,'DPR')
 if (tread==1) then
    anaddb_dtset%vs_qrad_tolkms(:) = dprarr(1:2)
    ABI_CHECK(dprarr(1) >= zero, "vs_qrad must be >= 0")
    ABI_CHECK(dprarr(2) > zero, "vs_tolkms must be > zero")
 end if
!W

!X

!Y

!Z

!=====================================================================
!end non-dependent variables
!=====================================================================

!=======================================================================
!Read in dependent variables (dependent on dimensions above)
!=======================================================================

!A

 ABI_ALLOCATE(anaddb_dtset%atifc,(natom))
 anaddb_dtset%atifc(:) = 0
 if(anaddb_dtset%natifc>=1)then
   ! default to 1 for first natifc atoms
   anaddb_dtset%atifc(1:anaddb_dtset%natifc)=1

   if(anaddb_dtset%natifc>marr)then
     marr=anaddb_dtset%natifc
     ABI_DEALLOCATE(intarr)
     ABI_DEALLOCATE(dprarr)
     ABI_ALLOCATE(intarr,(marr))
     ABI_ALLOCATE(dprarr,(marr))
   end if
   call intagm(dprarr,intarr,jdtset,marr,anaddb_dtset%natifc,string(1:lenstr),'atifc',tread,'INT')
   if(tread==1) anaddb_dtset%atifc(1:anaddb_dtset%natifc)= intarr(1:anaddb_dtset%natifc)
!  check of whether values of atifc are valid is done in chkin9
 end if

!B

!C

!D

!E

!F

!G

!H

!I

 ABI_ALLOCATE(anaddb_dtset%iatfix,(natom))
 anaddb_dtset%iatfix(:) = 0
 if ((anaddb_dtset%relaxat == 1).and.(anaddb_dtset%natfix > 0)) then
   if(natom > marr)then
     marr = natom
     ABI_DEALLOCATE(intarr)
     ABI_DEALLOCATE(dprarr)
     ABI_ALLOCATE(intarr,(marr))
     ABI_ALLOCATE(dprarr,(marr))
   end if
   call intagm(dprarr,intarr,jdtset,marr,anaddb_dtset%natfix,string(1:lenstr),'iatfix',tread,'INT')
   if(tread==1) anaddb_dtset%iatfix(1:anaddb_dtset%natfix) = intarr(1:anaddb_dtset%natfix)
 end if
!FIXME: need a test on values of iatfix: are they just 1 or 0?

 if ((anaddb_dtset%relaxstr == 1).and.(anaddb_dtset%nstrfix > 0)) then
   anaddb_dtset%istrfix(:) = 0
   call intagm(dprarr,intarr,jdtset,marr,anaddb_dtset%nstrfix,string(1:lenstr),'istrfix',tread,'INT')
   if(tread==1) anaddb_dtset%istrfix(1:anaddb_dtset%nstrfix) = intarr(1:anaddb_dtset%nstrfix)
 end if
!FIXME: need a test on values of istrfix

 if (anaddb_dtset%natprj_bs > 0) then
   ABI_ALLOCATE(anaddb_dtset%iatprj_bs,(anaddb_dtset%natprj_bs))
   if(anaddb_dtset%natprj_bs>marr)then
     marr=anaddb_dtset%natprj_bs
     ABI_DEALLOCATE(intarr)
     ABI_DEALLOCATE(dprarr)
     ABI_ALLOCATE(intarr,(marr))
     ABI_ALLOCATE(dprarr,(marr))
   end if
   anaddb_dtset%iatprj_bs(:)=0
   call intagm(dprarr,intarr,jdtset,marr,anaddb_dtset%natprj_bs,string(1:lenstr),'iatprj_bs',tread,'INT')
   if(tread==1) then
     anaddb_dtset%iatprj_bs(1:anaddb_dtset%natprj_bs)=intarr(1:anaddb_dtset%natprj_bs)
   else
     write(message,'(3a)')&
     'natprj_bs is non zero but iatprj_bs is absent ',ch10,&
     'Action: specify iatprj_bs in your input file.'
     MSG_ERROR(message)
   end if
 end if

!J

!K

!L

!M

!N

!O

!P

!Q

 if (anaddb_dtset%nqshft/=0)then
   if(3*anaddb_dtset%nqshft>marr)then
     marr=3*anaddb_dtset%nqshft
     ABI_DEALLOCATE(intarr)
     ABI_DEALLOCATE(dprarr)
     ABI_ALLOCATE(intarr,(marr))
     ABI_ALLOCATE(dprarr,(marr))
   end if
   anaddb_dtset%q1shft(:,:)=zero
   call intagm(dprarr,intarr,jdtset,marr,3*anaddb_dtset%nqshft, string(1:lenstr),'q1shft',tread,'DPR')
   if(tread==1) anaddb_dtset%q1shft(1:3,1:anaddb_dtset%nqshft)=&
&   reshape(dprarr(1:3*anaddb_dtset%nqshft),(/3,anaddb_dtset%nqshft/))
 end if

 ABI_ALLOCATE(anaddb_dtset%qph1l,(3,anaddb_dtset%nph1l))
 ABI_ALLOCATE(anaddb_dtset%qnrml1,(anaddb_dtset%nph1l))
 if (anaddb_dtset%nph1l/=0)then
   if(4*anaddb_dtset%nph1l>marr)then
     marr=4*anaddb_dtset%nph1l
     ABI_DEALLOCATE(intarr)
     ABI_DEALLOCATE(dprarr)
     ABI_ALLOCATE(intarr,(marr))
     ABI_ALLOCATE(dprarr,(marr))
   end if
   anaddb_dtset%qph1l(:,:)=zero
   anaddb_dtset%qnrml1(:)=zero
   call intagm(dprarr,intarr,jdtset,marr,4*anaddb_dtset%nph1l,string(1:lenstr),'qph1l',tread,'DPR')
   if(tread==1)then
     do iph1=1,anaddb_dtset%nph1l
       do ii=1,3
         anaddb_dtset%qph1l(ii,iph1)=dprarr(ii+(iph1-1)*4)
       end do
       anaddb_dtset%qnrml1(iph1)=dprarr(4+(iph1-1)*4)
       if(abs(anaddb_dtset%qnrml1(iph1))<DDB_QTOL)then
         write(message, '(5a)' )&
         'The first list of wavevectors ','should not have non-analytical data.',ch10,&
         'Action: correct the first list',' of wavevectors in the input file.'
         MSG_ERROR(message)
       end if
     end do
   end if
 end if

 ABI_ALLOCATE(anaddb_dtset%qph2l,(3,anaddb_dtset%nph2l))
 ABI_ALLOCATE(anaddb_dtset%qnrml2,(anaddb_dtset%nph2l))
 if (anaddb_dtset%nph2l/=0)then
   if(4*anaddb_dtset%nph2l>marr)then
     marr=4*anaddb_dtset%nph2l
     ABI_DEALLOCATE(intarr)
     ABI_DEALLOCATE(dprarr)
     ABI_ALLOCATE(intarr,(marr))
     ABI_ALLOCATE(dprarr,(marr))
   end if
   anaddb_dtset%qph2l(:,:)=zero
   anaddb_dtset%qnrml2(:)=zero
   call intagm(dprarr,intarr,jdtset,marr,4*anaddb_dtset%nph2l,string(1:lenstr),'qph2l',tread,'DPR')
   if(tread==1)then
     do iph2=1,anaddb_dtset%nph2l
       do ii=1,3
         anaddb_dtset%qph2l(ii,iph2)=dprarr(ii+(iph2-1)*4)
       end do
       anaddb_dtset%qnrml2(iph2)=dprarr(4+(iph2-1)*4)
       if(abs(anaddb_dtset%qnrml2(iph2))>DDB_QTOL)then
         write(message, '(5a)' )&
         'The second list of wavevectors',' should have only non-analytical data.',ch10,&
         'Action: correct the second list','of wavevectors in the input file.'
         MSG_ERROR(message)
       end if
     end do
   end if
 end if

 if (anaddb_dtset%nqpath > 0) then
   ABI_ALLOCATE(anaddb_dtset%qpath,(3,anaddb_dtset%nqpath))
   if(3*anaddb_dtset%nqpath>marr)then
     marr=3*anaddb_dtset%nqpath
     ABI_DEALLOCATE(intarr)
     ABI_DEALLOCATE(dprarr)
     ABI_ALLOCATE(intarr,(marr))
     ABI_ALLOCATE(dprarr,(marr))
   end if
   anaddb_dtset%qpath(:,:)=zero
   call intagm(dprarr,intarr,jdtset,marr,3*anaddb_dtset%nqpath,string(1:lenstr),'qpath',tread,'DPR')
   if(tread==1) then
     anaddb_dtset%qpath(1:3,1:anaddb_dtset%nqpath)= reshape(dprarr(1:3*anaddb_dtset%nqpath),(/3,anaddb_dtset%nqpath/))
   else
     write(message,'(3a)')&
     'nqpath is non zero but qpath is absent ',ch10,&
     'Action: specify qpath in your input file.'
     MSG_ERROR(message)
   end if
 end if

!R

!S

!T

!U

!V

!W

!X

!Y

!Z

!=======================================================================
!Finished reading in variables - deallocate
!=======================================================================

 ABI_DEALLOCATE(dprarr)
 ABI_DEALLOCATE(intarr)

!=======================================================================
!Check consistency of input variables:
!=======================================================================

 if (anaddb_dtset%frmin > anaddb_dtset%frmax) then
   write(message,'(3a)' )&
   'frmax should be higher than frmin',ch10,&
   'Action: change frmax and/or frmin  in your input file.'
   MSG_ERROR(message)
 end if

 if (anaddb_dtset%nqpath==0 .and. anaddb_dtset%elphflag==1) then
   write(message,'(4a)' )&
   'elphflag is 1 but no nqpath has been specified','for phonon linewidths',ch10,&
   'Action: specify nqpath and qpath(3,nqpath) in your input file.'
   MSG_ERROR(message)
 end if

 if(anaddb_dtset%telphint /= 2 .and. (anaddb_dtset%ep_b_min /= 0 .or. anaddb_dtset%ep_b_max /= 0)) then
   write(message, '(a,i0,3a)' )&
   'telphint is ',anaddb_dtset%telphint,', but ep_b_min or ep_b_max',ch10,&
   'are set /= 1. They will not be used'
   call wrtout(std_out,message,'COLL')
   MSG_WARNING(message)

 else if(anaddb_dtset%telphint == 2 .and. (anaddb_dtset%ep_b_min == 0 .or. anaddb_dtset%ep_b_max == 0)) then
   write(message, '(a,i0,4a)' )&
   'telphint is ',anaddb_dtset%telphint,', but ep_b_min or ep_b_max',ch10,&
   'are not both set. ',ch10,&
   'Action: set ep_b_min and ep_b_max in your input file.',ch10
   MSG_ERROR(message)
 end if

 if(anaddb_dtset%thmflag < 3) then
   if ((anaddb_dtset%telphint == 0 .or. anaddb_dtset%prtnest == 1 .or. &
        anaddb_dtset%prtnest == 2 .or. anaddb_dtset%prtfsurf== 1) .and. sum(anaddb_dtset%kptrlatt) == 0 ) then
     write (message, '(3a)') &
     'if tetrahedron integration is used, ',&
     'or the output of the nesting function/Fermi surface is required, ',&
     'you must specify the kptrlatt'
     MSG_ERROR(message)
   end if
 end if

 if(anaddb_dtset%prtdos/=0 .and. anaddb_dtset%ifcflag/=1) then
   write(message, '(3a)' )&
   'ifcflag must be 1 when the calculation of the phonon DOS is required ',ch10,&
   'Action: correct ifcflag in your input file.'
   MSG_ERROR(message)
 end if

 if(anaddb_dtset%prtsrlr/=0 .and. anaddb_dtset%ifcflag/=1) then
   write(message, '(3a)' )&
   'ifcflag must be 1 for the SR/LR decomposition of the phonon frequencies',ch10,&
   'Action: correct ifcflag in your input file.'
   MSG_ERROR(message)
 end if

 if (anaddb_dtset%gruns_nddbs /=0 .and. anaddb_dtset%ifcflag /=1) then
   MSG_ERROR("ifcflag must be 1 for Grunesein calculation")
 end if

 if (anaddb_dtset%vs_qrad_tolkms(1) /= zero .and. anaddb_dtset%ifcflag /=1) then
   MSG_ERROR("ifcflag must be 1 to calculate speed of sound")
 end if

 if(anaddb_dtset%prtdos/=0 .and. sum(abs(anaddb_dtset%ng2qpt(:))) < 3 ) then
   write(message, '(3a)' )&
   'ng2qpt must be specified when the calculation of the phonon DOS is required ',ch10,&
   'Action: correct ng2qpt in your input file.'
   MSG_ERROR(message)
 end if

 if (anaddb_dtset%ifltransport /= 0 .and. anaddb_dtset%ep_keepbands /= 1) then
   write(message, '(3a)' )&
   'Band dependency of electron phonon matrix elements must be kept for transport ',ch10,&
   'Action: set ep_keepbands to 1 in your input file.'
   MSG_ERROR(message)
 end if

 if (anaddb_dtset%ifltransport > 1 .and. sum(abs(anaddb_dtset%kptrlatt)) == 0) then
   write(message, '(3a)' )&
   'For inelastic transport or electron lifetime calculations you must specify kprtlatt ',ch10,&
   'Action: copy kptrlatt from your abinit GS file to your anaddb input file.'
   MSG_ERROR(message)
 end if

!FIXME: add check that if freeze_displ /= 0 then you need to be doing ifc and phonon interpolation

 if (anaddb_dtset%ifcflag > 0 .and. sum(abs(anaddb_dtset%ngqpt)) == 0) then
   write(message, '(3a)' )&
   'if you want interatomic force constant output, anaddb needs ngqpt input variable ',ch10,&
   'Action: set ngqpt in your input file.'
   MSG_ERROR(message)
 end if

!check that q-grid refinement is a divisor of ngqpt in each direction
 if(any(anaddb_dtset%qrefine(1:3) > 1) .and. &
    any(abs(dmod(dble(anaddb_dtset%ngqpt(1:3))/dble(anaddb_dtset%qrefine(1:3)),one)) > tol10) ) then
   write(message, '(a,3i10,a,a,a,3i8,a,a)' )&
   'qrefine is',anaddb_dtset%qrefine(1:3),' The only allowed values',ch10,&
   'are integers which are divisors of the ngqpt grid', anaddb_dtset%ngqpt(1:3),ch10,&
   'Action: correct qrefine in your input file.'
   MSG_ERROR(message)
 end if

!check that fermie and nelect are not both specified
 if(abs(anaddb_dtset%elph_fermie) > tol10 .and. abs(anaddb_dtset%ep_extrael) > tol10) then
   write(message, '(a,E10.2,a,E10.2,a,a,a)' )&
    'elph_fermie (',anaddb_dtset%elph_fermie,') and ep_extrael (',anaddb_dtset%ep_extrael, '), may not both be non 0',ch10,&
    'Action: remove one of the two in your input file.'
   MSG_ERROR(message)
 end if

 ! Check for possible typos.
 call anaddb_chkvars(string)

end subroutine invars9
!!***

!----------------------------------------------------------------------

!!****f* m_anaddb_dataset/outvars_anaddb
!!
!! NAME
!! outvars_anaddb
!!
!! FUNCTION
!! Open input file for the anaddb code, then
!! echoes the input information.
!!
!! INPUTS
!! anaddb_dtset= (derived datatype) contains all the input variables
!! nunit=unit number for input or output
!!
!! OUTPUT
!!  (only writing)
!!
!! NOTES
!! Should be executed by one processor only.
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!      chkvars_in_string,inupper
!!
!! SOURCE

subroutine outvars_anaddb (anaddb_dtset,nunit)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: nunit
 type(anaddb_dataset_type),intent(in) :: anaddb_dtset

!Local variables -------------------------
!scalars
 integer :: ii,iph1,iph2,iqpt,iqshft

!*********************************************************************

!Write the heading
 write(nunit,'(a,80a,a)') ch10,('=',ii=1,80),ch10
 write(nunit, '(2a)' )' -outvars_anaddb: echo values of input variables ----------------------',ch10

!The flags
 if (anaddb_dtset%dieflag/=0 .or. anaddb_dtset%ifcflag/=0 .or. &
     anaddb_dtset%nlflag/=0 .or. anaddb_dtset%thmflag/=0 .or. &
     anaddb_dtset%elaflag/=0 .or. anaddb_dtset%elphflag/=0 .or. &
     anaddb_dtset%polflag/=0 .or. anaddb_dtset%instrflag/=0 .or. &
     anaddb_dtset%piezoflag/=0) then
   write(nunit,'(a)')' Flags :'
   if(anaddb_dtset%dieflag/=0)write(nunit,'(3x,a9,3i10)')'  dieflag',anaddb_dtset%dieflag
   if(anaddb_dtset%ifcflag/=0)write(nunit,'(3x,a9,3i10)')'  ifcflag',anaddb_dtset%ifcflag
   if(anaddb_dtset%nlflag/=0)write(nunit,'(3x,a9,3i10)')'   nlflag',anaddb_dtset%nlflag
   if(anaddb_dtset%thmflag/=0)write(nunit,'(3x,a9,3i10)')'  thmflag',anaddb_dtset%thmflag
   if(anaddb_dtset%elaflag/=0)write(nunit,'(3x,a9,3i10)')'  elaflag',anaddb_dtset%elaflag
   if(anaddb_dtset%elphflag/=0)write(nunit,'(3x,a9,3i10)')' elphflag',anaddb_dtset%elphflag
   if(anaddb_dtset%polflag/=0)write(nunit,'(3x,a9,3i10)')'  polflag',anaddb_dtset%polflag
   if(anaddb_dtset%instrflag/=0)write(nunit,'(3x,a9,3i10)')'instrflag',anaddb_dtset%instrflag
   if(anaddb_dtset%piezoflag/=0)write(nunit,'(3x,a9,3i10)')'piezoflag',anaddb_dtset%piezoflag
 end if

!Write the general information
 if (anaddb_dtset%rfmeth/=1 .or. &
     anaddb_dtset%enunit/=0 .or. &
     anaddb_dtset%eivec/=0 .or. &
     anaddb_dtset%asr/=0 .or. &
     anaddb_dtset%chneut/=0 .or. &
     anaddb_dtset%selectz/=0 .or. anaddb_dtset%symdynmat/=1) then
   write(nunit,'(a)')' Miscellaneous information :'
   if(anaddb_dtset%rfmeth/=1)write(nunit,'(3x,a9,3i10)')'   rfmeth',anaddb_dtset%rfmeth
   if(anaddb_dtset%enunit/=0)write(nunit,'(3x,a9,3i10)')'   enunit',anaddb_dtset%enunit
   if(anaddb_dtset%eivec/=0) write(nunit,'(3x,a9,3i10)')'    eivec',anaddb_dtset%eivec
   if(anaddb_dtset%asr/=0)   write(nunit,'(3x,a9,3i10)')'      asr',anaddb_dtset%asr
   if(anaddb_dtset%chneut/=0)write(nunit,'(3x,a9,3i10)')'   chneut',anaddb_dtset%chneut
   if(anaddb_dtset%selectz/=0)write(nunit,'(3x,a9,3i10)')'  selectz',anaddb_dtset%selectz
   if(anaddb_dtset%symdynmat/=1)write(nunit,'(3x,a9,3i10)')'symdynmat',anaddb_dtset%symdynmat
 end if
 if(anaddb_dtset%prtvol/=0) write(nunit,'(3x,a9,i10)')'   prtvol',anaddb_dtset%prtvol

!Frequency information
 if(anaddb_dtset%dieflag==1)then
   write(nunit,'(a)')' Frequency information :'
   write(nunit,'(3x,a9,3i10)')'    nfreq',anaddb_dtset%nfreq
   write(nunit,'(3x,a9,7x,3es16.8)')'    frmin',anaddb_dtset%frmin
   write(nunit,'(3x,a9,7x,3es16.8)')'    frmax',anaddb_dtset%frmax
 end if

!For interatomic force constant information
 if(anaddb_dtset%ifcflag/=0)then
   write(nunit,'(a)')' Interatomic Force Constants Inputs :'
   write(nunit,'(3x,a9,3i10)')'   dipdip',anaddb_dtset%dipdip
   if(anaddb_dtset%nsphere/=0)write(nunit,'(3x,a9,3i10)')'  nsphere',anaddb_dtset%nsphere
   if(abs(anaddb_dtset%rifcsph)>tol10)write(nunit,'(3x,a9,E16.6)')'  nsphere',anaddb_dtset%rifcsph
   write(nunit,'(3x,a9,3i10)')'   ifcana',anaddb_dtset%ifcana
   write(nunit,'(3x,a9,3i10)')'   ifcout',anaddb_dtset%ifcout
   if(anaddb_dtset%natifc>=1)then
     write(nunit,'(3x,a9,3i10)')'   natifc',anaddb_dtset%natifc
     write(nunit,'(3x,a9,8i10)')'    atifc',(anaddb_dtset%atifc(ii),ii=1,anaddb_dtset%natifc)
   end if
   write(nunit,'(a)')' Description of grid 1 :'
   write(nunit,'(3x,a9,3i10)')'     brav',anaddb_dtset%brav
   write(nunit,'(3x,a9,3i10)')'    ngqpt',anaddb_dtset%ngqpt(1:3)
   write(nunit,'(3x,a9,3i10)')'   nqshft',anaddb_dtset%nqshft
   if (anaddb_dtset%nqshft/=0)then
     write(nunit,'(3x,a9)')'   q1shft'
     do iqshft=1,anaddb_dtset%nqshft
       write(nunit,'(19x,4es16.8)') (anaddb_dtset%q1shft(ii,iqshft),ii=1,3)
     end do
   end if
   if (any(anaddb_dtset%qrefine(:) > 1)) then
     write(nunit,'(3x,a9,3i10)')'  qrefine', anaddb_dtset%qrefine
   end if
   ! Speed of sound
   if (anaddb_dtset%vs_qrad_tolkms(1) > zero) then
      write(nunit,'(a,2es16.8)')"vs_qrad_tolkms", (anaddb_dtset%vs_qrad_tolkms(:))
   end if
 end if

!Phonon density of states with gaussian method
 if(anaddb_dtset%prtdos/=0)then
   write(nunit,'(a)')' Phonon DOS information :'
   write(nunit,'(3x,a9,es16.8)')'dosdeltae',anaddb_dtset%dosdeltae
   write(nunit,'(3x,a9,es16.8)')' dossmear',anaddb_dtset%dossmear
 end if

!Thermal information
 if(anaddb_dtset%thmflag/=0)then
   write(nunit,'(a)')' Thermal information :'
   write(nunit,'(3x,a9,3i10)')'    nchan',anaddb_dtset%nchan
   write(nunit,'(3x,a9,3i10)')'   nwchan',anaddb_dtset%nwchan
   write(nunit,'(3x,a9,7x,3es16.8)')'   dostol',anaddb_dtset%dostol
   write(nunit,'(3x,a9,7x,3es16.8)')'   thmtol',anaddb_dtset%thmtol
   write(nunit,'(3x,a9,3i10)')'  ntemper',anaddb_dtset%ntemper
   write(nunit,'(3x,a9,7x,3es16.8)')'temperinc',anaddb_dtset%temperinc
   write(nunit,'(3x,a9,7x,3es16.8)')'tempermin',anaddb_dtset%tempermin
 endif

!Grid 2 description
 if(anaddb_dtset%thmflag/=0 .or. anaddb_dtset%prtdos/=0)then
   write(nunit,'(a)')' Description of grid 2 (Fourier interp. or BZ sampling):'
   write(nunit,'(3x,a9,3i10)')'   ng2qpt',anaddb_dtset%ng2qpt(1:3)
   write(nunit,'(3x,a9,3i10)')'   ngrids',anaddb_dtset%ngrids
   write(nunit,'(3x,a9,7x,3es16.8)')'   q2shft',anaddb_dtset%q2shft(1:3)
 end if

!Non-linear response information
 if (anaddb_dtset%nlflag /= 0) then
   write(nunit,'(a)')' Non-linear response information :'
   write(nunit,'(3x,a9,i10)') '   alphon',anaddb_dtset%alphon
   write(nunit,'(3x,a9,3i10)')'   prtmbm',anaddb_dtset%prtmbm
   write(nunit,'(3x,a9,3i10)')'  ramansr',anaddb_dtset%ramansr
 end if

!Structural relaxation at fixed polarization
 if (anaddb_dtset%polflag /= 0) then
   write(nunit,'(a)')' Relaxation at fixed polarization :'
   if (anaddb_dtset%relaxat == 1) then
     write(nunit,'(3x,a9,i10)') '  relaxat',anaddb_dtset%relaxat
   end if
   if (anaddb_dtset%relaxstr == 1) then
     write(nunit,'(a12,i10)') ' relaxstr',anaddb_dtset%relaxstr
   end if
 end if

!Elphon information
 if (anaddb_dtset%elphflag /= 0) then
   write(nunit,'(a)')' Elphon calculation will be carried out'
   write(nunit,'(a12,E16.6)') 'elphsmear', anaddb_dtset%elphsmear
   write(nunit,'(a12,E16.6)') 'a2fsmear', anaddb_dtset%a2fsmear
   write(nunit,'(a12,E16.6)') 'mustar', anaddb_dtset%mustar
   write(nunit,'(a12,i10)') 'nqpath', anaddb_dtset%nqpath
   write(nunit,'(a12)') 'qpath'
   do iqpt=1,anaddb_dtset%nqpath
     write(nunit,'(12x,3(E16.6,1x))') anaddb_dtset%qpath(:,iqpt)
   end do
   write(nunit,'(a12,i10)') 'telphint', anaddb_dtset%telphint
   if (anaddb_dtset%telphint == 0) then
     write(nunit,'(a)') ' Tetrahedron integration for elphon'
   else if (anaddb_dtset%telphint == 1) then
     write(nunit,'(a)') ' Smeared weight integration for elphon'
   else if (anaddb_dtset%telphint == 2) then
     write(nunit,'(a)') ' Band filtered integration for elphon'
   end if
   if (abs(anaddb_dtset%elph_fermie) > tol10) then
     write(nunit,'(a12,E16.6)')  'elph_fermie', anaddb_dtset%elph_fermie
   end if
   if (anaddb_dtset%ep_extrael /= 0) then
     if (abs(anaddb_dtset%ep_extrael) > 1.0d2) then
        write(nunit,'(a,E20.12)')' Doping set by the user is (negative for el doping) :',anaddb_dtset%ep_extrael
     else
       write(nunit,'(a,E16.6)')  'Elphon: extra electrons per unit cell = ', anaddb_dtset%ep_extrael
     end if
   end if
   if (anaddb_dtset%ep_nspline /= 20) then
     write(nunit,'(a,I8)')  'Elphon: scale factor for spline interpolation in RTA = ', anaddb_dtset%ep_nspline
   end if
   if (anaddb_dtset%band_gap < 10.0d0) then
     write(nunit,'(a,E16.6)')  'Elphon: set band gap to (in eV) = ', anaddb_dtset%band_gap
   end if

   if (sum(abs(anaddb_dtset%kptrlatt)) > 0) then
     write(nunit,'(a12,3(3(i3,1x),2x))' ) 'kptrlatt',reshape( anaddb_dtset%kptrlatt(:,:), (/9/) )
   end if

   if (sum(abs(anaddb_dtset%kptrlatt_fine)) > 0) then
     write(nunit,'(a12,3(3(i3,1x),2x))' ) 'kptrlatt_fine ',reshape( anaddb_dtset%kptrlatt_fine(:,:), (/9/) )
   end if

   if (anaddb_dtset%ep_keepbands == 1) then
     write(nunit, '(a)') ' Will keep band dependency in gkk in memory.'
     write(nunit, '(a)') ' WARNING: the memory requirements will be multiplied by nbands**2 !!!'
   end if

   if (anaddb_dtset%ep_scalprod == 1) then
     write(nunit, '(a)') ' scalar product will be performed when assembling the gamma matrices.'
     write(nunit, '(a)') ' WARNING: with this option you can not distinguish which '
     write(nunit, '(a)') '    linewidth comes from which phonon mode !!!'
   end if

   if (anaddb_dtset%prtbltztrp== 1) write(nunit, '(a)') ' Will output input files for BoltzTraP'
   if (anaddb_dtset%prtfsurf == 1) write(nunit, '(a)') ' Will output fermi surface in XCrysDen format'
   if (anaddb_dtset%prt_ifc == 1) write(nunit, '(a)') ' Will output real space IFC in AI2PS and TDEP format'
   if (anaddb_dtset%prtnest == 1) write(nunit, '(a)') ' Will output nesting factor'

   if (anaddb_dtset%ifltransport == 1) then
     write(nunit, '(a)') ' Will perform transport calculation in elphon to get'
     write(nunit, '(a,a)') ' resistivity and thermal conductivity as a function of T',ch10
     write(nunit, '(a,es16.6,a)' ) ' Minimum temperature for transport outputs: ', anaddb_dtset%tempermin, ' K'
     write(nunit, '(a,es16.6,a)' ) ' Maximum temperature for transport outputs: ', &
       anaddb_dtset%tempermin+anaddb_dtset%temperinc*anaddb_dtset%ntemper, ' K'
     write(nunit, '(a,i6)' ) ' Number of temperature points for transport outputs: ', anaddb_dtset%ntemper
     write(nunit, '(a)' )
   end if

   if (anaddb_dtset%gkqwrite == 1) then
     write(nunit,'(a,a)' ) 'Gkk matrix elements on input grid of ',&
     'qpoints will be written to disk. File gkqfile must be absent.'
   end if
   if (anaddb_dtset%gkk_rptwrite == 1) then
     write(nunit,'(a,a)' ) 'Gkk matrix elements in real space ',&
     'will be written to disk. File gkk_rpt_file must be absent.'
   end if
   if (anaddb_dtset%gkk2write == 1) then
     write(nunit,'(a,a)' ) 'Full grid gkk matrix elements ',&
     'will be written to disk. File gkk2file must be absent.'
   end if
 end if

 if (anaddb_dtset%gruns_nddbs /= 0) then
   write(nunit,'(a)' ) "Will compute Gruneisen parameters with finite difference method. DDB files:"
   do ii=1,anaddb_dtset%gruns_nddbs
     write(nunit, "(2a)")"    ",trim(anaddb_dtset%gruns_ddbs(ii))
   end do
 end if

!List of vector 1  (reduced coordinates)
 if(anaddb_dtset%nph1l/=0)then
   write(nunit,'(a)')' First list of wavevector (reduced coord.) :'
   write(nunit,'(3x,a9,3i10)')'    nph1l',anaddb_dtset%nph1l
   write(nunit,'(3x,a9)')'    qph1l'
   do iph1=1,anaddb_dtset%nph1l
     write(nunit,'(19x,3es16.8,2x,es11.3)') &
       (anaddb_dtset%qph1l(ii,iph1),ii=1,3),anaddb_dtset%qnrml1(iph1)
   end do
 end if

!List of vector 2  (cartesian coordinates)
 if(anaddb_dtset%nph2l/=0)then
   write(nunit,'(a)')' Second list of wavevector (cart. coord.) :'
   write(nunit,'(3x,a9,3i10)')'    nph2l',anaddb_dtset%nph2l
   write(nunit,'(3x,a9)')'    qph2l'
   do iph2=1,anaddb_dtset%nph2l
     write(nunit,'(19x,3es16.8,2x,es11.3)') (anaddb_dtset%qph2l(ii,iph2),ii=1,3),anaddb_dtset%qnrml2(iph2)
   end do
 end if

!phonon frozen in supercell
 if (abs(anaddb_dtset%freeze_displ) > tol10) then
   write(nunit,'(a)') 'Phonon displacements will be output, frozen into supercells'
   write(nunit,'(a,E20.10)') ' Chosen amplitude of frozen displacements = ', anaddb_dtset%freeze_displ
 end if

!atom projected bs files
 if (abs(anaddb_dtset%natprj_bs) > 0) then
   write(nunit,'(a)') 'Phonon band structure files, with atomic projections, will be output '
   write(nunit,'(a)') ' Chosen atoms for projection = '
   write(nunit,'(10I6)') anaddb_dtset%iatprj_bs
 end if

 write(nunit,'(a,80a,a)') ch10,('=',ii=1,80),ch10

end subroutine outvars_anaddb
!!***

!----------------------------------------------------------------------

!!****f* m_anaddb_dataset/anaddb_init
!!
!! NAME
!! anaddb_init
!!
!! FUNCTION
!! Initialize the code ppddb9: write heading and make the first i/os
!!
!! INPUTS
!!  input_path: String with input file path. Empty string activates files file legacy mode.
!!
!! OUTPUT
!! character(len=fnlen) filnam(7)=character strings giving file names
!!
!! NOTES
!! 1. Should be executed by one processor only.
!! 2. File names refer to following files, in order:
!!     (1) Formatted input file
!!     (2) Formatted output file
!!     (3) Input Derivative Database
!!     (4) Output Molecular Dynamics
!!     (5) Input electron-phonon matrix elements
!!     (6) Root name for electron-phonon file names
!!     (7) Name of file containing the 3 ddk filenames and the GS wf file name
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!      chkvars_in_string,inupper
!!
!! SOURCE

subroutine anaddb_init(input_path, filnam)

!Arguments -------------------------------
!arrays
 character(len=*),intent(in) :: input_path
 character(len=*),intent(out) :: filnam(7)

!Local variables -------------------------
!scalars
 integer :: lenstr, marr, jdtset, tread
 character(len=strlen) :: string
!arrays
 integer,allocatable :: intarr(:)
 real(dp),allocatable :: dprarr(:)

! *********************************************************************

 if (len_trim(input_path) == 0) then
   ! Legacy Files file mode.
   write(std_out, "(2a)")" DeprecationWarning: ",ch10
   write(std_out, "(a)") "     The files file has been deprecated in Abinit9 and will be removed in Abinit10."
   write(std_out, "(2a)")"     Use the syntax `anaddb t01.abi` where t01.abi is an anaddb input with ddb_filepath.",ch10
   write(std_out, "(3a)")'            ddb_filepath = "out_DDB"',ch10,ch10

   write(std_out,*)' Give name for formatted input file: '
   read(std_in, '(a)' ) filnam(1)
   write(std_out,'(a,a)' )'-   ',trim(filnam(1))
   write(std_out,*)' Give name for formatted output file: '
   read(std_in, '(a)' ) filnam(2)
   write(std_out,'(a,a)' )'-   ',trim(filnam(2))
   write(std_out,*)' Give name for input derivative database: '
   read(std_in, '(a)' ) filnam(3)
   write(std_out,'(a,a)' )'-   ',trim(filnam(3))
   write(std_out,*)' Give name for output molecular dynamics: '
   read(std_in, '(a)' ) filnam(4)
   write(std_out,'(a,a)' )'-   ',trim(filnam(4))
   write(std_out,*)' Give name for input elphon matrix elements (GKK file): '
   read(std_in, '(a)' ) filnam(5)
   write(std_out,'(a,a)' )'-   ',trim(filnam(5))
   write(std_out,*)' Give root name for elphon output files: '
   read(std_in, '(a)' ) filnam(6)
   write(std_out,'(a,a)' )'-   ',trim(filnam(6))
   write(std_out,*)' Give name for file containing ddk filenames for elphon/transport: '
   read(std_in, '(a)' ) filnam(7)
   write(std_out,'(a,a)' )'-   ',trim(filnam(7))

 else
   ! Read input
   call instrng(input_path, lenstr, 1, strlen, string)
   ! To make case-insensitive, map characters to upper case.
   call inupper(string(1:lenstr))

   filnam = ""
   filnam(1) = input_path
   filnam(2) = "run.abo"

   marr = 3
   ABI_MALLOC(intarr, (marr))
   ABI_MALLOC(dprarr, (marr))
   jdtset = 0

   ! Allow user to override default values
   call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), "output_file", tread, 'KEY', key_value=filnam(2))
   write(std_out, "(2a)")'- Name for formatted output file: ', trim(filnam(2))

   call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), "ddb_filepath", tread, 'KEY', key_value=filnam(3))
   ABI_CHECK(tread == 1, "ddb_filepath variable must be specified in the input file")
   write(std_out, "(2a)")'- Input derivative database: ', trim(filnam(3))

   ! Nobody knows the scope of this line in the files file.
   !call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), "md_output", tread, 'KEY', key_value=filnam(4))
   call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), "gkk_filepath", tread, 'KEY', key_value=filnam(5))
   if (tread == 1) write(std_out, "(2a)")'- Name for input elphon matrix elements (GKK file): ', trim(filnam(5))

   call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), "eph_prefix", tread, 'KEY', key_value=filnam(6))
   if (tread == 1) write(std_out, "(2a)")"- Root name for elphon output files: ", trim(filnam(6))

   call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), "ddk_filepath", tread, 'KEY', key_value=filnam(7))
   if (tread == 1) write(std_out, "(2a)")"- File containing ddk filenames for elphon/transport: ", trim(filnam(7))

   ABI_FREE(intarr)
   ABI_FREE(dprarr)
 end if

end subroutine anaddb_init
!!***

!!****f* m_anaddb_dataset/anaddb_chkvars
!! NAME
!! anaddb_chkvars
!!
!! FUNCTION
!!  Examines the input string, to check whether all names are allowed.
!!
!! INPUTS
!!  string*(*)=string of character
!!   the string (with upper case) from the input file, to which the XYZ data is (possibly) appended
!!
!! OUTPUT
!!
!! PARENTS
!!      m_anaddb_dataset
!!
!! CHILDREN
!!      chkvars_in_string,inupper
!!
!! SOURCE

subroutine anaddb_chkvars(string)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: string

!Local variables-------------------------------
!scalars
 integer,parameter :: protocol0=0
 character(len=100) :: list_logicals,list_strings
 character(len=10000) :: list_vars

!************************************************************************

!Here, list all admitted variable names (max 10 per line, to fix the ideas)
!Note: Do not use "double quotation mark" for the string since it triggers a bug in docchk.py (abirules script)
!<ANADDB_VARS>
!A
 list_vars=                 ' alphon asr a2fsmear atifc'
!B
 list_vars=trim(list_vars)//' brav band_gap'
!C
 list_vars=trim(list_vars)//' chneut'
!D
 list_vars=trim(list_vars)//' dieflag dipdip dossum dosdeltae dossmear dostol'
!E
 list_vars=trim(list_vars)//' ep_scalprod eivec elaflag elphflag enunit'
 list_vars=trim(list_vars)//' ep_b_min ep_b_max ep_int_gkk ep_keepbands ep_nqpt ep_nspline ep_prt_yambo'
 list_vars=trim(list_vars)//' elphsmear elph_fermie ep_extrael ep_qptlist'
!F
 list_vars=trim(list_vars)//' freeze_displ frmax frmin'
!G
 list_vars=trim(list_vars)//' gkk2write gkk_rptwrite gkqwrite gruns_nddbs'
!H
!I
 list_vars=trim(list_vars)//' ifcana ifcflag ifcout ifltransport instrflag istrfix iatfix iatprj_bs'
!J
!K
 list_vars=trim(list_vars)//' kptrlatt kptrlatt_fine'
!L
!M
 list_vars=trim(list_vars)//' mustar'
!N
 list_vars=trim(list_vars)//' natfix natifc natom natprj_bs nchan ndivsm nfreq ngrids nlflag nph1l nph2l'
 list_vars=trim(list_vars)//' nqpath nqshft nsphere nstrfix ntemper nwchan ngqpt ng2qpt'
!O
 list_vars=trim(list_vars)//' outboltztrap'
!P
 list_vars=trim(list_vars)//' piezoflag polflag prtddb prtdos prt_ifc prtmbm prtfsurf'
 list_vars=trim(list_vars)//' prtnest prtphbands prtsrlr prtvol prtbltztrp'
!Q
 list_vars=trim(list_vars)//' qrefine qgrid_type q1shft q2shft qnrml1 qnrml2 qpath qph1l qph2l'
!R
 list_vars=trim(list_vars)//' ramansr relaxat relaxstr rfmeth rifcsph'
!S
 list_vars=trim(list_vars)//' selectz symdynmat symgkq'
!T
 list_vars=trim(list_vars)//' targetpol telphint thmflag temperinc tempermin thermal_supercell thmtol'
!U
 list_vars=trim(list_vars)//' use_k_fine'
!V
 list_vars=trim(list_vars)//' vs_qrad_tolkms'
!W
!X
!Y
!Z

!Logical input variables
 list_logicals=' '

!String input variables
 list_strings=' gruns_ddbs ddb_filepath output_file gkk_filepath eph_prefix ddk_filepath' ! md_output
!</ANADDB_VARS>

!Extra token, also admitted:
!<ANADDB_UNITS>
 list_vars=trim(list_vars)//' au Angstr Angstrom Angstroms Bohr Bohrs eV Ha'
 list_vars=trim(list_vars)//' Hartree Hartrees K nm Ry Rydberg Rydbergs T Tesla'
!</ANADDB_UNITS>

!<ANADDB_OPERATORS>
 list_vars=trim(list_vars)//' sqrt '
!</ANADDB_OPERATORS>

!Transform to upper case
 call inupper(list_vars)
 call inupper(list_logicals)
 call inupper(list_strings)

 call chkvars_in_string(protocol0, list_vars, list_logicals, list_strings, string)

end subroutine anaddb_chkvars
!!***

end module m_anaddb_dataset
!!***

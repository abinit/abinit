!!****m*ABINIT/m_anaddb_dataset
!! NAME
!!  m_anaddb_dataset
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2025 ABINIT group (XG,JCC,CL,MVeithen,XW,MJV)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_anaddb_dataset

 use defs_basis
 use m_xmpi
 use m_abicore
 use m_errors
 use m_nctk
 use netcdf

 use m_parser,    only : instrng
 use m_fstrings,  only : next_token, rmquotes, sjoin, inupper, ltoa, itoa, basename
 use m_clib,      only : clib_mkdir_if_needed
 use m_matrix,    only : mati3det
 use m_parser,    only : intagm, chkvars_in_string, instrng
 use m_crystal,   only : crystal_t
 use m_ddb,       only : DDB_QTOL, chkin9
 use m_ddb_hdr,   only : ddb_hdr_type

 implicit none

 private

 public:: anaddb_dataset_type
!!***

!----------------------------------------------------------------------

!!****t*m_anaddb_dataset/anaddb_dataset_type
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
  integer:: alphon
  integer:: asr
  integer:: brav
  integer:: chneut
  integer:: dieflag
  integer:: dipdip
  integer:: dipquad
  integer:: dossum
  integer:: dos_maxmode
  integer:: ep_scalprod
  integer:: eivec
  integer:: elaflag
  integer:: elphflag
  integer:: enunit
  integer:: flexoflag
  integer:: gkk2write
  integer:: gkk_rptwrite
  integer:: gkqwrite
  integer:: gruns_nddbs
  integer:: ifcana
  integer:: ifcflag
  integer:: ifcout
  integer:: ifltransport
  integer:: instrflag
  integer:: lwf_anchor_proj
  integer:: lwf_disentangle
  integer:: lwf_nwann
  integer:: lwfflag
  integer:: natfix
  integer:: natifc
  integer:: natprj_bs
  integer:: nchan
  integer:: ndivsm = 20
  integer:: nfreq
  integer:: ngrids
  integer:: nlflag
  integer:: nph1l
  integer:: nph2l
  integer:: nqpath
  integer:: nqshft
  integer:: nsphere
  integer:: nstrfix
  integer:: ntemper
  integer:: nwchan
  integer:: outboltztrap
  integer:: piezoflag
  integer:: polflag
  integer:: prtdos
  integer:: prt_ifc
  integer:: prtddb
  integer:: prtmbm
  integer:: prtfsurf
  integer:: prtnest
  integer:: prtphbands
  integer:: prtsrlr  ! print the short-range/long-range decomposition of phonon freq.
  integer:: prtvol = 0
  integer:: ramansr
  integer:: relaxat
  integer:: relaxstr
  integer:: rfmeth
  integer:: selectz
  integer:: symdynmat
  integer:: telphint
  integer:: thmflag
  integer:: qgrid_type
  integer:: quadquad
  integer:: ep_b_min
  integer:: ep_b_max
  integer:: ep_int_gkk
  integer:: ep_keepbands
  integer:: ep_nqpt
  integer:: ep_nspline
  integer:: ep_prt_yambo
  integer:: symgkq
  integer:: use_k_fine
  integer:: prtbltztrp

  ! These are not input variables, but important dimensions read from ddb
  integer:: natom
  integer:: msize
  integer:: mpert
  integer:: ntypat
  integer:: usepaw  ! GA: TODO Remove usepaw

  integer:: ngqpt(9)             ! ngqpt(9) instead of ngqpt(3) is needed in wght9.f
  integer:: istrfix(6)
  integer:: lwf_ngqpt(3)
  integer:: ng2qpt(3)
  integer:: qrefine(3)
  integer:: kptrlatt(3, 3)
  integer:: kptrlatt_fine(3, 3)
  integer:: thermal_supercell(3, 3)

! Real(dp)
  real(dp):: a2fsmear
  real(dp):: band_gap
  real(dp):: dosdeltae
  real(dp):: dossmear
  real(dp):: dostol
  real(dp):: elphsmear
  real(dp):: elph_fermie
  real(dp):: ep_extrael
  real(dp):: freeze_displ
  real(dp):: frmax
  real(dp):: frmin

  real(dp):: lwf_anchor_qpt(3)
  real(dp):: lwf_mu
  real(dp):: lwf_sigma
  real(dp):: temperinc
  real(dp):: tempermin
  real(dp):: thmtol
  real(dp):: mustar
  real(dp):: rifcsph

  real(dp):: q1shft(3, 4)
  real(dp):: q2shft(3)
  real(dp):: targetpol(3)
  real(dp):: vs_qrad_tolkms(2) = 0

  character(len=fnlen):: filename_input
  character(len=fnlen):: filename_output
  character(len=fnlen):: prefix_outdata
  character(len=fnlen):: filename_ddb
  character(len=fnlen):: filename_ddk
  character(len=fnlen):: prefix_eph
  character(len=fnlen):: filename_gkk
  character(len=fnlen):: filename_eigr2d

  character(len=strlen):: input_string
  ! The entire input string.
  integer:: lenstr  ! Length of the entire input string.

! Integer arrays
  integer, allocatable:: atifc(:)
   ! atifc(natom) Atoms for which IFC should be computed.

  integer, allocatable:: atifcflg(:)
   ! atifcflg(natom) Flag (0 or 1) to tell which atom is analysed in ifc.

  integer, allocatable:: iatfix(:)
  ! iatfix(natom)

  integer, allocatable:: iatprj_bs(:)

  integer, allocatable:: lwf_anchor_iband(:)
  integer, allocatable:: lwf_projector(:)

! Real arrays
  real(dp), allocatable:: qnrml1(:)
  ! qnrml1(nph1l)

  real(dp), allocatable:: qnrml2(:)
  ! qnrml2(nph2l)

  real(dp), allocatable:: qpath(:,:)
  ! qpath(3, nqpath)

  real(dp), allocatable:: qph1l(:,:)
  ! qph1l(3, nph1l)

  real(dp), allocatable:: qph2l(:,:)
  ! qph2l(3, nph2l)

  real(dp), allocatable:: ep_qptlist(:,:)
  ! qph2l(3, ep_nqpt)

  character(len = fnlen), allocatable:: gruns_ddbs(:)
  ! gruns_ddbs(gruns_nddbs)

  contains

    procedure :: init => anaddb_dtset_init
     ! Construct the object from the dtset.

    procedure :: free => anaddb_dtset_free
     ! Free dynamic memory.

    procedure :: read_input => anaddb_dtset_read_input
     ! Read input file, and some info from the ddb or IFC.

    procedure :: bcast_files => anaddb_dtset_bcast_files
     ! Broadcast file names

    procedure :: outvars => outvars_anaddb
     ! Broadcast file names

 end type anaddb_dataset_type
!!***

contains
!!***

!!****f*m_anaddb_dataset/anaddb_dtset_free
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
!! SOURCE

subroutine anaddb_dtset_free(dtset)

!Arguments------------------------------------
!scalars
 class(anaddb_dataset_type), intent(inout):: dtset

! *************************************************************************

 ABI_SFREE(dtset%atifc)
 ABI_SFREE(dtset%atifcflg)
 ABI_SFREE(dtset%iatfix)
 ABI_SFREE(dtset%iatprj_bs)
 ABI_SFREE(dtset%qnrml1)
 ABI_SFREE(dtset%qnrml2)
 ABI_SFREE(dtset%qpath)
 ABI_SFREE(dtset%qph1l)
 ABI_SFREE(dtset%qph2l)
 ABI_SFREE(dtset%ep_qptlist)
 ABI_SFREE(dtset%gruns_ddbs)
 if (dtset%lwfflag==1) then
    ABI_SFREE(dtset%lwf_anchor_iband)
 else if (dtset%lwfflag==2) then
    ABI_SFREE(dtset%lwf_projector)
 end if

end subroutine anaddb_dtset_free
!!***

!----------------------------------------------------------------------

!!****f*m_anaddb_dataset/invars9
!!
!! NAME
!! invars9
!!
!! FUNCTION
!! Open input file for the anaddb code, then reads or echoes the input information.
!!
!! INPUTS
!! lenstr = actual length of string
!! natom = number of atoms, needed for atifc
!! string*(*)=string of characters containing all input variables and data
!!
!! OUTPUT
!! dtset= (derived datatype) contains all the input variables
!!
!! NOTES
!! Should be executed by one processor only.
!!
!! 27/01/2009: MJV: I have cleaned this routine extensively, putting all
!!  variables in alphabetical order, and in a second segment the dependent
!!  variables which need to be allocated depending on the dimensions read in.
!!  Could be divided into two routines as in abinit.
!!  FIXME: move checks to chkin9?
!!
!! SOURCE

subroutine invars9(dtset, lenstr, natom, string)

!Arguments-------------------------------
!scalars
 integer, intent(in):: lenstr, natom
 character(len=*), intent(in):: string
 type(anaddb_dataset_type), intent(inout):: dtset

!Local variables-------------------------
!scalars
 integer, parameter:: vrsddb = 100401  ! Set routine version number here:
 integer, parameter:: jdtset = 1
 integer:: ii, iph1, iph2, marr, tread, start, idet
 character(len = 500):: message
 character(len = fnlen):: path
!arrays
 integer, allocatable:: intarr(:)
 real(dp), allocatable:: dprarr(:)

!*********************************************************************
 marr = 3
 ABI_MALLOC(intarr, (marr))
 ABI_MALLOC(dprarr, (marr))

!copy natom to dtset
 dtset%natom = natom

!=====================================================================
!start reading in dimensions and non-dependent variables
!=====================================================================

!A

!typical value for gaussian smearing of a2F function
 dtset%a2fsmear = 0.00002_dp
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'a2fsmear',tread, 'ENE')
 if(tread == 1) dtset%a2fsmear = dprarr(1)
 if (dtset%a2fsmear < tol6) then
   write(message, '(a, f10.3, a, a, a, a, a)' )&
   'a2fsmear is ',dtset%a2fsmear, ', but only values > 1.e-6 ',ch10, &
   'are allowed',ch10, 'Action: correct a2fsmear in your input file.'
   ABI_ERROR(message)
 end if

 dtset%alphon = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'alphon',tread, 'INT')
 if(tread == 1) dtset%alphon = intarr(1)
!FIXME: need a test on input value

 dtset%asr = 1
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'asr',tread, 'INT')
 if(tread == 1) dtset%asr = intarr(1)
 if(dtset%asr < -2 .or. dtset%asr > 5)then
   write(message, '(a, i0, 5a)' )&
   'asr is ',dtset%asr, ', but the only allowed values',ch10, &
   'are 0, 1, 2, 3, 4, 5, -1 or-2 .',ch10, 'Action: correct asr in your input file.'
!  Note : negative values are allowed when the acoustic sum rule
!  is to be applied after the analysis of IFCs
!  3, 4 are for rotational invariance (under development)
!  5 is for hermitian imposition of the ASR
   ABI_ERROR(message)
 end if

!B

!Target band gap in eV
!The default value is just a very large number that will not be used in changing the band gap
 dtset%band_gap = 999.0d0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'band_gap',tread, 'DPR')
 if(tread == 1) dtset%band_gap = dprarr(1)

 dtset%brav = 1
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'brav',tread, 'INT')
 if(tread == 1) dtset%brav = intarr(1)
 if(dtset%brav <= -2 .or. dtset%brav >= 5 .or. dtset%brav == 0)then
   write(message, '(a, i0, a5)' )&
   'brav is ',dtset%brav, ', but the only allowed values',ch10, &
   'are-1, 1, 2, 3 or 4 .',ch10, 'Action: correct brav in your input file.'
   ABI_ERROR(message)
 end if

!C

 dtset%chneut = 1
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'chneut',tread, 'INT')
 if(tread == 1) dtset%chneut = intarr(1)
 if(dtset%chneut < 0 .or. dtset%chneut > 2)then
   write(message, '(a, i0, 5a)' )&
   'chneut is ',dtset%chneut, ', but the only allowed values',ch10, &
   'are 0, 1 or 2.',ch10, 'Action: correct chneut in your input file.'
   ABI_ERROR(message)
 end if

!D

 dtset%dieflag = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'dieflag',tread, 'INT')
 if(tread == 1) dtset%dieflag = intarr(1)
 if(dtset%dieflag < 0 .or. dtset%dieflag > 4)then
   write(message, '(a, i0, 5a)' )&
   'dieflag is ',dtset%dieflag, ', but the only allowed values',ch10, &
   'are 0, 1, 2, 3 or 4.',ch10, 'Action: correct dieflag in your input file.'
   ABI_ERROR(message)
 end if

 dtset%dipdip = 1
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'dipdip',tread, 'INT')
 if(tread == 1) dtset%dipdip = intarr(1)
 if(dtset%dipdip < -1 .or. dtset%dipdip > 1)then
   write(message, '(a, i0, 5a)' )&
   'dipdip is ',dtset%dipdip, ', but the only allowed values',ch10, &
   'are-1, 0 or 1 .',ch10, 'Action: correct dipdip in your input file.'
   ABI_ERROR(message)
 end if

 dtset%dipquad = 1
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'dipquad',tread, 'INT')
 if(tread == 1) dtset%dipquad = intarr(1)
 if(dtset%dipquad < -1 .or. dtset%dipquad > 1)then
   write(message, '(a, i0, 5a)' )&
   'dipquad is ',dtset%dipquad, ', but the only allowed values',ch10, &
   'are 0 or 1 .',ch10, 'Action: correct dipquad in your input file.'
   ABI_ERROR(message)
 end if

 dtset%ep_scalprod = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'ep_scalprod',tread, 'INT')
 if(tread == 1) dtset%ep_scalprod = intarr(1)
 if(dtset%ep_scalprod < 0 .or. dtset%ep_scalprod > 1) then
   write(message, '(a, i0, 5a)' )&
   'ep_scalprod is ',dtset%ep_scalprod, ', but the only allowed values',ch10, &
   'are 0 or 1.',ch10, 'Action: correct ep_scalprod in your input file.'
   ABI_ERROR(message)
 end if

 dtset%dosdeltae = 0.2_dp/Ha_cmm1
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'dosdeltae',tread, 'DPR')
 if(tread == 1) dtset%dosdeltae = dprarr(1)

!FIXME : should probably be smaller
 dtset%dossmear = one/Ha_cmm1
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'dossmear',tread, 'DPR')
 if(tread == 1) dtset%dossmear = dprarr(1)
 if(dtset%dossmear <= zero)then
   write(message, '(a, es14.4, 3a)' )&
   'dossmear is ',dtset%dossmear, ', which is lower than 0 .',ch10, &
   'Action: correct dossmear in your input file.'
   ABI_ERROR(message)
 end if

 dtset%dostol = 0.25_dp
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'dostol',tread, 'DPR')
 if(tread == 1) dtset%dostol = dprarr(1)
 if(dtset%dostol < zero)then
   write(message, '(a, es14.4, 3a)' )&
   'dostol is ',dtset%dostol, ', which is lower than 0 .',ch10, &
   'Action: correct dostol in your input file.'
   ABI_ERROR(message)
 end if

 dtset%dossum = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'dossum',tread, 'INT')
 if(tread == 1) dtset%dossum = intarr(1)
 if(dtset%dossum < 0 .or. dtset%dossum > one)then
   write(message, '(a, i0, 5a)' )&
   'dossum is ',dtset%dossum, ', but the only allowed values',ch10, &
   'are 0, 1',ch10, 'Action: correct dossum in your input file.'
   ABI_ERROR(message)
 end if

 dtset%dos_maxmode = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'dos_maxmode',tread, 'INT')
 if(tread == 1) dtset%dos_maxmode = intarr(1)
 if(dtset%dos_maxmode < 0 .or. dtset%dos_maxmode > 3*natom)then
   write(message, '(a, i0, 5a)' )&
   'dos_maxmode is ',dtset%dos_maxmode, ', but the only allowed values',ch10, &
   'are 0 to 3*natom ',ch10, 'Action: correct dos_maxmode in your input file.'
   ABI_ERROR(message)
 end if

!E

 dtset%eivec = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'eivec',tread, 'INT')
 if(tread == 1) dtset%eivec = intarr(1)
 if(dtset%eivec < 0 .or. dtset%eivec > 4)then
   write(message, '(a, i0, 5a)' )&
   'eivec is ',dtset%eivec, ', but the only allowed values',ch10, &
   'are 0, 1, 2, 3 or 4.',ch10, 'Action: correct eivec in your input file.'
   ABI_ERROR(message)
 end if

 dtset%elaflag = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'elaflag',tread, 'INT')
 if(tread == 1) dtset%elaflag = intarr(1)
 if(dtset%elaflag < 0 .or. dtset%elaflag > 5)then
   write(message, '(a, i0, 5a)' )&
   'elaflag is ',dtset%elaflag, ', but the only allowed values',ch10, &
   'are 0, 1, 2, 3, 4 or 5 .',ch10, 'Action: correct elaflag in your input file.'
   ABI_ERROR(message)
 end if

!By default use the real fermie (tests for abs(elph_fermie) < tol10 in the code)
 dtset%elph_fermie = zero
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'elph_fermie',tread, 'ENE')
 if(tread == 1) dtset%elph_fermie = dprarr(1)

!extra charge in unit cell (number of electrons) wrt neutral cell
!holes are negative values (reduce number of electrons)
 dtset%ep_extrael = zero
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'ep_extrael',tread, 'DPR')
 if(tread == 1) dtset%ep_extrael = dprarr(1)

!number to control the spline interpolation in RTA
 dtset%ep_nspline = 20
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'ep_nspline',tread, 'INT')
 if(tread == 1) dtset%ep_nspline = intarr(1)
 if(dtset%ep_nspline < 0 .or. dtset%ep_nspline > 1000) then
   write(message, '(a, i0, 5a)' )&
   'ep_nspline is ',dtset%ep_nspline, ', but this should not be ',ch10, &
   'negative or too large .',ch10, 'Action: correct ep_nspline in your input file.'
   ABI_ERROR(message)
 end if

!interpolate gkk or gamma. It should be better to interpolate gkk onto the
!k_phon, since the integration weights will be treated the same way
 dtset%ep_int_gkk = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'ep_int_gkk',tread, 'INT')
 if(tread == 1) dtset%ep_int_gkk = intarr(1)
 if(dtset%ep_int_gkk < 0 .or. dtset%ep_int_gkk > 1) then
   write(message, '(a, i0, 5a)' )&
   'ep_int_gkk is ',dtset%ep_int_gkk, ', but the only allowed values',ch10, &
   'are 0 or 1.',ch10, 'Action: correct ep_int_gkk in your input file.'
   ABI_ERROR(message)
 end if

 dtset%elphflag = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'elphflag',tread, 'INT')
 if(tread == 1) dtset%elphflag = intarr(1)
 if(dtset%elphflag < 0 .or. dtset%elphflag > 1)then
   write(message, '(a, i0, 5a)' )&
   'elphflag = ',dtset%elphflag, ', but the allowed values',ch10, &
   'are 0, or 1.',ch10, 'Action: correct elphflag in your input file.'
   ABI_ERROR(message)
 end if

!typical value for gaussian smearing, but can vary sensibly with the metal
 dtset%elphsmear = 0.01_dp
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'elphsmear',tread, 'ENE')
 if(tread == 1) dtset%elphsmear = dprarr(1)
 if (dtset%elphsmear < tol6) then
   write(message, '(a, f10.3, 5a)' )&
   'elphsmear is ',dtset%elphsmear, '. Only values > 1.e-6 ',ch10, &
   'are allowed',ch10, 'Action: correct elphsmear in your input file.'
   ABI_ERROR(message)
 end if

 dtset%enunit = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'enunit',tread, 'INT')
 if(tread == 1) dtset%enunit = intarr(1)
 if(dtset%enunit < 0 .or. dtset%enunit > 2)then
   write(message, '(a, i0, 5a)' )&
   'enunit is ',dtset%enunit, ', but the only allowed values',ch10, &
   'are 0, 1 or 2.',ch10, 'Action: correct enunit in your input file.'
   ABI_ERROR(message)
 end if

!Default is 0-not used unless telphint == 2
 dtset%ep_b_max = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'ep_b_max',tread, 'INT')
 if(tread == 1) then
   dtset%ep_b_max = intarr(1)
   if(dtset%ep_b_max < 1) then
     write(message, '(a, i0, 5a)' )&
     'ep_b_max is ',dtset%ep_b_max, ', but the only allowed values',ch10, &
     'are between 1 and nband.',ch10, 'Action: correct ep_b_max in your input file.'
     ABI_ERROR(message)
   end if
 end if

!Default is 0-not used unless telphint == 2
 dtset%ep_b_min = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'ep_b_min',tread, 'INT')
 if(tread == 1) then
   dtset%ep_b_min = intarr(1)
   if(dtset%ep_b_min < 1) then
     write(message, '(a, i0, 5a)' )&
     'ep_b_min is ',dtset%ep_b_min, ', but the only allowed values',ch10, &
     'are between 1 and nband.',ch10, 'Action: correct ep_b_min in your input file.'
     ABI_ERROR(message)
   end if
 end if

 dtset%ep_keepbands = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'ep_keepbands',tread, 'INT')
 if(tread == 1) dtset%ep_keepbands = intarr(1)
 if(dtset%ep_keepbands < 0 .or. dtset%ep_keepbands > 1) then
   write(message, '(a, i0, 5a)' )&
   'ep_keepbands is ',dtset%ep_keepbands, ', but the only allowed values',ch10, &
   'are 0 or 1 .',ch10, 'Action: correct ep_keepbands in your input file.'
   ABI_ERROR(message)
 end if

 dtset%ep_nqpt = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'ep_nqpt',tread, 'INT')
 if(tread == 1) dtset%ep_nqpt = intarr(1)
 if(dtset%ep_nqpt < 0) then
   write(message, '(a, i0, 5a)' )&
   'ep_nqpt is ',dtset%ep_nqpt, ', but the only allowed values',ch10, &
   'are > 0.',ch10, 'Action: correct ep_nqpt in your input file.'
   ABI_ERROR(message)
 end if

 if (dtset%ep_nqpt > 0) then
   ABI_MALLOC(dtset%ep_qptlist, (3, dtset%ep_nqpt))
   if(3*dtset%ep_nqpt > marr)then
     marr = 3*dtset%ep_nqpt
     ABI_FREE(intarr)
     ABI_FREE(dprarr)
     ABI_MALLOC(intarr, (marr))
     ABI_MALLOC(dprarr, (marr))
   end if
   dtset%ep_qptlist(:,:)=zero
   call intagm(dprarr, intarr, jdtset, marr, 3*dtset%ep_nqpt, string(1:lenstr), 'ep_qptlist',tread, 'DPR')
   if(tread == 1) then
     dtset%ep_qptlist(1:3, 1:dtset%ep_nqpt)=&
     reshape(dprarr(1:3*dtset%ep_nqpt), (/3, dtset%ep_nqpt/))
   else
     write(message, '(3a)')&
     'ep_nqpt is non zero but ep_qptlist is absent ',ch10, &
     'Action: specify ep_qptlist in your input file.'
     ABI_ERROR(message)
   end if
 end if


!F

 dtset%flexoflag = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'flexoflag',tread, 'INT')
 if(tread == 1) dtset%flexoflag = intarr(1)
 if(dtset%flexoflag < 0 .or. dtset%flexoflag > 4)then
   write(message, '(3a, i0, 5a)' )&
   ' flexoflag is ',dtset%flexoflag, ', but the only allowed values',ch10, &
   'are 0, 1, 2, 3, 4  .',ch10, 'Action: correct flexoflag in your input file.'
   ABI_ERROR(message)
 end if

 dtset%freeze_displ = zero
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'freeze_displ',tread, 'DPR')
 if(tread == 1) dtset%freeze_displ = dprarr(1)


 dtset%frmax = ten
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'frmax',tread, 'DPR')
 if(tread == 1) dtset%frmax = dprarr(1)
 if (dtset%frmax < 0) then
   write(message, '(a, f10.3, 5a)' )&
   'frmax is ',dtset%frmax, '. Only values > 0 ',ch10, &
   'are allowed',ch10, 'Action: correct frmax in your input file.'
   ABI_ERROR(message)
 end if

 dtset%frmin = zero
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'frmin',tread, 'DPR')
 if(tread == 1) dtset%frmin = dprarr(1)
 if (dtset%frmin < 0) then
   write(message, '(a, f10.3, 5a)' )&
   'frmin is ',dtset%frmin, '. Only values > 0 ',ch10, &
   'are allowed',ch10, 'Action: correct frmin in your input file.'
   ABI_ERROR(message)
 end if

!G

 dtset%gkk2write = 0
!call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'gkk2write',&
!& tread, 'INT')
!if(tread == 1) dtset%gkk2write = intarr(1)
!if(dtset%gkk2write < 0 .or. dtset%gkk2write > 1) then
!write(message, '(a, a, a, i8, a, a, a, a, a)' )&
!&  ' invars9 : ERROR -',ch10, &
!&  '  gkk2write is',dtset%gkk2write, &
!&  ', but the only allowed values',ch10, &
!&  '  are 0 or 1 .',ch10, &
!&  '  Action: correct gkk2write in your input file.'
!ABI_ERROR(message)
!end if

 dtset%gkk_rptwrite = 0
!call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'gkk_rptwrite',&
!& tread, 'INT')
!if(tread == 1) dtset%gkk_rptwrite = intarr(1)
!if(dtset%gkk_rptwrite < 0 .or. dtset%gkk_rptwrite > 1) then
!write(message, '(a, a, a, i8, a, a, a, a, a)' )&
!&  ' invars9 : ERROR -',ch10, &
!&  '  gkk_rptwrite is',dtset%gkk_rptwrite, &
!&  ', but the only allowed values',ch10, &
!&  '  are 0 or 1 .',ch10, &
!&  '  Action: correct gkk_rptwrite in your input file.'
!ABI_ERROR(message)
!end if

 dtset%gkqwrite = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'gkqwrite',tread, 'INT')
 if(tread == 1) dtset%gkqwrite = intarr(1)
 if(dtset%gkqwrite < 0 .or. dtset%gkqwrite > 1) then
   write(message, '(a, i0, 5a)' )&
   'gkqwrite is ',dtset%gkqwrite, ', but the only allowed values',ch10, &
   'are 0 or 1.',ch10, 'Action: correct gkqwrite in your input file.'
   ABI_ERROR(message)
 end if

 dtset%gruns_nddbs = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'gruns_nddbs',tread, 'INT')
 if (tread == 1) dtset%gruns_nddbs = intarr(1)

 if (dtset%gruns_nddbs /= 0) then
   ! Read list of DDB paths.
   ABI_MALLOC(dtset%gruns_ddbs, (dtset%gruns_nddbs))
   start = index(string, "GRUNS_DDBS") + len("GRUNS_DDBS") + 1
   do ii = 1, dtset%gruns_nddbs
     if (next_token(string, start, path) /= 0) then
       ABI_ERROR(sjoin("Cannot find DDB path in input string:", ch10, string(start:)))
     end if
     dtset%gruns_ddbs(ii) = rmquotes(path)
   end do
 end if

!H

!I

 dtset%ifcana = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'ifcana',tread, 'INT')
 if(tread == 1) dtset%ifcana = intarr(1)
 if(dtset%ifcana < 0 .or. dtset%ifcana > 1)then
   write(message, '(a, i0, 5a)' )&
   'ifcana is ',dtset%ifcana, ', but the only allowed values',ch10, &
   'are 0 or 1.',ch10, 'Action: correct ifcana in your input file.'
   ABI_ERROR(message)
 end if

 dtset%ifcflag = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'ifcflag',tread, 'INT')
 if(tread == 1) dtset%ifcflag = intarr(1)
 if(dtset%ifcflag < 0 .or. dtset%ifcflag > 1)then
   write(message, '(a, i0, 5a)' )&
   'ifcflag is ',dtset%ifcflag, ', but the only allowed values',ch10, &
   'are 0 or 1.',ch10, 'Action: correct ifcflag in your input file.'
   ABI_ERROR(message)
 end if

 dtset%ifcout = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'ifcout',tread, 'INT')
 if(tread == 1) dtset%ifcout = intarr(1)
 if(dtset%ifcout < -1)then
   write(message, '(a, i0, 3a)' )&
   'ifcout is ',dtset%ifcout, ', which is lower than-1.',ch10, &
   'Action: correct ifcout in your input file.'
   ABI_ERROR(message)
 end if

 dtset%ifltransport = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'ifltransport',tread, 'INT')
 if(tread == 1) dtset%ifltransport = intarr(1)
 if(dtset%ifltransport < 0 .or. dtset%ifltransport > 3) then
   write(message, '(a, i0, 5a)' )&
   'ifltransport is ',dtset%ifltransport, ', but the only allowed values',ch10, &
   'are 0 or 1 or 2 or 3.',ch10, 'Action: correct ifltransport in your input file.'
   ABI_ERROR(message)
 end if

 dtset%instrflag = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'instrflag',tread, 'INT')
 if(tread == 1) dtset%instrflag = intarr(1)
 if(dtset%instrflag < 0 .or. dtset%instrflag > 1)then
   write(message, '(a, i0, 5a)' )&
   'instrflag is ',dtset%instrflag, ', but the only allowed values',ch10, &
   'are 0, 1.',ch10, 'Action: correct instrflag in your input file.'
   ABI_ERROR(message)
 end if

!J

!K

 dtset%kptrlatt = 0
!why this test on reading in kptrlatt?
 marr = 9
 ABI_FREE(intarr)
 ABI_FREE(dprarr)
 ABI_MALLOC(intarr, (marr))
 ABI_MALLOC(dprarr, (marr))
 call intagm(dprarr, intarr, jdtset, marr, 9, string(1:lenstr), 'kptrlatt',tread, 'INT')
 if(tread == 1)dtset%kptrlatt(1:3, 1:3)=reshape(intarr(1:9), (/3, 3/))
!NOTE: no a priori way to test the validity of the integers in kptrlatt

 dtset%kptrlatt_fine(:,:)=0
 marr = 9
 ABI_FREE(intarr)
 ABI_FREE(dprarr)
 ABI_MALLOC(intarr, (marr))
 ABI_MALLOC(dprarr, (marr))
 call intagm(dprarr, intarr, jdtset, marr, 9, string(1:lenstr), 'kptrlatt_fine',tread, 'INT')
 if(tread == 1)dtset%kptrlatt_fine(1:3, 1:3)=reshape(intarr(1:9), (/3, 3/))


!L

! lwf_*: lattic wannier function
dtset%lwf_anchor_qpt(:) = 0.0_dp
call intagm(dprarr, intarr, jdtset, marr, 3, string(1:lenstr), 'lwf_anchor_qpt',tread, 'DPR')
if(tread == 1) dtset%lwf_anchor_qpt(1:3) = dprarr(1:3)

dtset%lwf_disentangle = 0
call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'lwf_disentangle',tread, 'INT')
if(tread == 1) dtset%lwf_disentangle = intarr(1)
if(dtset%lwf_disentangle < 0 .or. dtset%lwf_disentangle > 3)then
   write(message, '(a, i0, 5a)' )&
        'lwf_disentangle is ',dtset%lwf_disentangle, ', but the only allowed values',ch10, &
        'are 0, 1.',ch10, 'Action: correct lwf_disentangle in your input file.'
   ABI_ERROR(message)
end if

dtset%lwf_anchor_proj = 0
call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'lwf_anchor_proj',tread, 'INT')
if(tread == 1) dtset%lwf_anchor_proj = intarr(1)
if(dtset%lwf_anchor_proj < 0 .or. dtset%lwf_anchor_proj > 2)then
   write(message, '(a, i0, 5a)' )&
        'lwf_anchor_proj is ',dtset%lwf_anchor_proj, ', but the only allowed values',ch10, &
        'are 0, 1, 2.',ch10, 'Action: correct lwf_anchor_proj in your input file.'
   ABI_ERROR(message)
end if


dtset%lwf_ngqpt(:)=0
call intagm(dprarr, intarr, jdtset, marr, 3, string(1:lenstr), 'lwf_ngqpt',tread, 'INT')
if(tread == 1) dtset%lwf_ngqpt(1:3)=intarr(1:3)
do ii = 1, 3
   if(dtset%lwf_ngqpt(ii)<0)then
      write(message, '(a, i0, a, i0, 3a, i0, a)' )&
           'lwf_ngqpt(',ii, ') is ',dtset%lwf_ngqpt(ii), ', which is lower than 0 .',ch10, &
           'Action: correct lwf_ngqpt(',ii, ') in your input file.'
      ABI_ERROR(message)
   end if
end do


dtset%lwfflag = 0
call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'lwfflag',tread, 'INT')
if(tread == 1) dtset%lwfflag = intarr(1)
if(dtset%lwfflag < 0 .or. dtset%lwfflag > 2)then
   write(message, '(a, i0, 5a)' )&
        'lwfflag is ',dtset%lwfflag, ', but the only allowed values',ch10, &
        'are 0, 1 and 2.',ch10, 'Action: correct lwfflag in your input file.'
   ABI_ERROR(message)
end if


dtset%lwf_mu = 0.0_dp
call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'lwf_mu',tread, 'DPR')
if(tread == 1) dtset%lwf_mu = dprarr(1)


dtset%lwf_nwann = 0
call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'lwf_nwann',tread, 'INT')
if(tread == 1) dtset%lwf_nwann = intarr(1)
if(dtset%lwf_nwann > natom*3)then
   write(message, '(a, i0, 2a, i0, 3a)' )&
        'lwf_nwann is ',dtset%lwf_nwann, ', which is larger than natom*3',' (=',natom*3, ')',ch10, &
        'Action: correct lwf_nwann in your input file.'
   ABI_ERROR(message)
end if

if( dtset%lwfflag > 0 .and. dtset%lwf_nwann .le. 0)then
   write(message, '(a, i0, 3a)' )&
        'lwf_nwann is ',dtset%lwf_nwann, ', which is not positive',ch10, &
        'Action: correct lwf_nwann in your input file.'
   ABI_ERROR(message)
end if

dtset%lwf_sigma = 0.01_dp
call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'lwf_sigma',tread, 'DPR')
if(tread == 1) dtset%lwf_sigma = dprarr(1)

!M

!typical value for mustar, but can vary sensibly with the metal
 dtset%mustar = 0.1_dp
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'mustar',tread, 'DPR')
 if(tread == 1) dtset%mustar = dprarr(1)
 if (dtset%mustar < zero) then
   write(message, '(a, f10.3, 5a)' )&
   'mustar is ',dtset%mustar, ', but only positive values',ch10, &
   'are allowed',ch10, 'Action: correct mustar in your input file.'
   ABI_ERROR(message)
 end if

!N

 dtset%natfix = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'natfix',tread, 'INT')
 if(tread == 1) dtset%natfix = intarr(1)
 if(dtset%natfix > natom)then
   write(message, '(a, i0, 2a, i0, 3a)' )&
   'natfix is ',dtset%natfix, ', which is larger than natom',' (=',natom, ')',ch10, &
   'Action: correct natfix in your input file.'
   ABI_ERROR(message)
 end if

 if(dtset%natfix < 0)then
   write(message, '(a, i0, 3a)' )&
   'natfix is ',dtset%natfix, ', which is < 0',ch10, &
   'Action: correct natfix in your input file.'
   ABI_ERROR(message)
 end if

 dtset%natifc = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'natifc',tread, 'INT')
 if(tread == 1) dtset%natifc = intarr(1)
 if(dtset%natifc < 0)then
   write(message, '(a, i0, 3a)' )&
   'natifc is ',dtset%natifc, ', which is lower than 0 .',ch10, &
   'Action: correct natifc in your input file.'
   ABI_ERROR(message)
 end if

 dtset%natprj_bs = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'natprj_bs',tread, 'INT')
 if(tread == 1) dtset%natprj_bs = intarr(1)
 if(dtset%natprj_bs < 0 .or. dtset%natprj_bs > natom)then
   write(message, '(a, i0, a, i0, 2a)' )&
   'natprj_bs is ',dtset%natprj_bs, ', but must be between 0 and natom = ',natom, ch10, &
   'Action: correct natprj_bs in your input file.'
   ABI_ERROR(message)
 end if

 dtset%nchan = 800
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'nchan',tread, 'INT')
 if(tread == 1) dtset%nchan = intarr(1)
!FIXME: check this-it should probably be .ge. 1, not 0
 if(dtset%nchan < 0)then
   write(message, '(a, i0, 3a)' )&
   'nchan is ',dtset%nchan, ', which is lower than 0 .',ch10, &
   'Action: correct nchan in your input file.'
   ABI_ERROR(message)
 end if

 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'ndivsm',tread, 'INT')
 if(tread == 1) dtset%ndivsm = intarr(1)
 if(dtset%ndivsm <= 0)then
   write(message, '(a, i0, 3a)' )&
   'ndivsm is ',dtset%ndivsm, ', which is <= 0 .',ch10, &
   'Action: correct ndivsm in your input file.'
   ABI_ERROR(message)
 end if

 dtset%nfreq = 1
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'nfreq',tread, 'INT')
 if(tread == 1) dtset%nfreq = intarr(1)
 if(dtset%nfreq < 0)then
   write(message, '(a, i0, 3a)' )&
   'nfreq is ',dtset%nfreq, ', which is lower than 0 .',ch10, &
   'Action: correct nfreq in your input file.'
   ABI_ERROR(message)
 end if

 dtset%ng2qpt(:)=0
 call intagm(dprarr, intarr, jdtset, marr, 3, string(1:lenstr), 'ng2qpt',tread, 'INT')
 if(tread == 1) dtset%ng2qpt(:)=intarr(1:3)
 do ii = 1, 3
   if(dtset%ng2qpt(ii)<0)then
     write(message, '(a, i0, a, i0, 3a, i0, a)' )&
     'ng2qpt(',ii, ') is ',dtset%ng2qpt(ii), ', which is lower than 0 .',ch10, &
     'Action: correct ng2qpt(',ii, ') in your input file.'
     ABI_ERROR(message)
   end if
 end do

 dtset%ngqpt(:)=0
 call intagm(dprarr, intarr, jdtset, marr, 3, string(1:lenstr), 'ngqpt',tread, 'INT')
 if(tread == 1) dtset%ngqpt(1:3)=intarr(1:3)
 do ii = 1, 3
   if(dtset%ngqpt(ii)<0)then
     write(message, '(a, i0, a, i0, 3a, i0, a)' )&
     'ngqpt(',ii, ') is ',dtset%ngqpt(ii), ', which is lower than 0 .',ch10, &
     'Action: correct ngqpt(',ii, ') in your input file.'
     ABI_ERROR(message)
   end if
 end do

 dtset%ngrids = 4
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'ngrids',tread, 'INT')
 if(tread == 1) dtset%ngrids = intarr(1)
 if(dtset%ngrids < 0)then
   write(message, '(a, i0, 3a)' )&
   'ngrids is ',dtset%ngrids, ', which is lower than 0 .',ch10, &
   'Action: correct ngrids in your input file.'
   ABI_ERROR(message)
 end if

 dtset%nlflag = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'nlflag',tread, 'INT')
 if(tread == 1) dtset%nlflag = intarr(1)
 if(dtset%nlflag < 0 .or. dtset%nlflag > 3)then
   write(message, '(a, i0, 5a)' )&
   'nlflag is ',dtset%nlflag, ', but the only allowed values',ch10, &
   'are 0, 1, 2 or 3.',ch10, 'Action: correct nlflag in your input file.'
   ABI_ERROR(message)
 end if

 dtset%nph1l = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'nph1l',tread, 'INT')
 if(tread == 1) dtset%nph1l = intarr(1)
 if(dtset%nph1l < 0)then
   write(message, '(a, i0, 3a)' )&
   'nph1l is ',dtset%nph1l, ', which is lower than 0 .',ch10, &
   'Action: correct nph1l in your input file.'
   ABI_ERROR(message)
 end if

 dtset%nph2l = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'nph2l',tread, 'INT')
 if(tread == 1) dtset%nph2l = intarr(1)
 if(dtset%nph2l < 0)then
   write(message, '(a, i0, 3a)' )&
   'nph2l is ',dtset%nph2l, ', which is lower than 0 .',ch10, &
   'Action: correct nph2l in your input file.'
   ABI_ERROR(message)
 end if

 dtset%nqpath = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'nqpath',tread, 'INT')
 if(tread == 1) dtset%nqpath = intarr(1)
 if(dtset%nqpath < 0)then
   write(message, '(a, i0, 3a)' )&
   'nqpath is ',dtset%nqpath, ', but must be positive',ch10, &
   'Action: correct elphflag in your input file.'
   ABI_ERROR(message)
 end if

 dtset%nqshft = 1
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'nqshft',tread, 'INT')
 if(tread == 1) dtset%nqshft = intarr(1)
 if(dtset%nqshft < 0 .or. dtset%nqshft == 3 .or. dtset%nqshft >= 5 )then
   write(message, '(a, i0, 5a)' )&
   'nqshft is ',dtset%nqshft, ', but the only allowed values',ch10, &
   'are 1, 2 or 4 .',ch10, 'Action: correct nqshft in your input file.'
   ABI_ERROR(message)
 end if

 dtset%nsphere = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'nsphere',tread, 'INT')
 if(tread == 1) dtset%nsphere = intarr(1)
 if(dtset%nsphere < -1)then
   write(message, '(a, i0, 3a)' )&
   'nsphere is ',dtset%nsphere, ', while it must be >= 0 or equal to-1',ch10, &
   'Action: correct nsphere in your input file.'
   ABI_ERROR(message)
 end if

 dtset%nstrfix = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'nstrfix',tread, 'INT')
 if(tread == 1) dtset%nstrfix = intarr(1)
 if(dtset%nstrfix > 6)then
   write(message, '(a, i0, 3a)' )&
   'nstrfix is ',dtset%nstrfix, ', which is larger than 6',ch10, &
   'Action: correct nstrfix in your input file.'
   ABI_ERROR(message)
 end if

 if(dtset%nstrfix < 0)then
   write(message, '(a, i0, 3a)' )&
   'nstrfix is ',dtset%nstrfix, ', which is < 0',ch10, &
   'Action: correct nstrfix in your input file.'
   ABI_ERROR(message)
 end if

 dtset%ntemper = 10
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'ntemper',tread, 'INT')
 if(tread == 1) dtset%ntemper = intarr(1)
 if(dtset%ntemper < 0)then
   write(message, '(a, i0, 3a)' )&
   'ntemper is ',dtset%ntemper, ', which is lower than 0',ch10, &
   'Action: correct ntemper in your input file.'
   ABI_ERROR(message)
 end if

 dtset%nwchan = 10
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'nwchan',tread, 'INT')
 if(tread == 1) dtset%nwchan = intarr(1)
!FIXME: check this-it should probably be .ge. 1, not 0
 if(dtset%nwchan < 0)then
   write(message, '(a, i0, 3a)' )&
   'nwchan is ',dtset%nwchan, ', which is lower than 0 .',ch10, &
   'Action: correct nwchan in your input file.'
   ABI_ERROR(message)
 end if

!O
 dtset%outboltztrap = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'outboltztrap',tread, 'INT')
 if(tread == 1) dtset%outboltztrap = intarr(1)
 if(dtset%outboltztrap < 0 .or. dtset%outboltztrap > 1)then
   write(message, '(a, i0, 5a)' )&
   'outboltztrap is ',dtset%outboltztrap, ', but the only allowed values',ch10, &
   'are 0 or 1.',ch10, 'Action: correct outboltztrap in your input file.'
   ABI_ERROR(message)
 end if


!P
 dtset%piezoflag = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'piezoflag',tread, 'INT')
 if(tread == 1) dtset%piezoflag = intarr(1)
 if(dtset%piezoflag < 0 .or. dtset%piezoflag > 7)then
   write(message, '(3a, i0, 5a)' )&
   ' piezoflag is ',dtset%piezoflag, ', but the only allowed values',ch10, &
   'are 0, 1, 2, 3, 4, 5, 6, 7  .',ch10, 'Action: correct piezoflag in your input file.'
   ABI_ERROR(message)
 end if

 dtset%polflag = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'polflag',tread, 'INT')
 if(tread == 1) dtset%polflag = intarr(1)
 if(dtset%polflag < 0 .or. dtset%polflag > 1)then
   write(message, '(a, i0, 5a)' )&
   'polflag is ',dtset%polflag, ', but the only allowed values',ch10, &
   'are 0 or 1.',ch10, 'Action: correct polflag in your input file.'
   ABI_ERROR(message)
 end if

 dtset%prtbltztrp = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'prtbltztrp',tread, 'INT')
 if(tread == 1) dtset%prtbltztrp = intarr(1)
 if(dtset%prtbltztrp < 0 .or. dtset%prtbltztrp > 1)then
   write(message, '(a, i0, 5a)' )&
   'prtbltztrp is ',dtset%prtbltztrp, ', but the only allowed values',ch10, &
   'are 0 or 1.',ch10, 'Action: correct prtbltztrp in your input file.'
   ABI_ERROR(message)
 end if

 ! Default is no output for PHDOS
 dtset%prtdos = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'prtdos',tread, 'INT')
 if(tread == 1) dtset%prtdos = intarr(1)
 if(dtset%prtdos < 0 .or. dtset%prtdos > 2) then
   write(message, '(a, i0, 5a)' )&
   'prtdos is ',dtset%prtdos, ', but the only allowed values',ch10, &
   'are 0 (no output) or 1 (gaussians) or 2 (tetrahedra) ',ch10, &
   'Action: correct prtdos in your input file.'
   ABI_ERROR(message)
 end if

!Default is no output for the Fermi Surface
 dtset%prtfsurf = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'prtfsurf',tread, 'INT')
 if(tread == 1) dtset%prtfsurf = intarr(1)
 if(dtset%prtfsurf < 0 .or. dtset%prtfsurf > 2) then
   write(message, '(a, i0, 5a)' )&
   'prtfsurf is ',dtset%prtfsurf, '. The only allowed values',ch10, &
   'are 0 (no output) or 1 (Xcrysden bxsf format)',ch10,  &
   'Action: correct prtfsurf in your input file.'
   ABI_ERROR(message)
 end if

!Default is no output of the real space IFC to file
 dtset%prt_ifc = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'prt_ifc',tread, 'INT')
 if(tread == 1) dtset%prt_ifc = intarr(1)
 if(dtset%prt_ifc < 0 .or. dtset%prt_ifc > 1) then
   write(message, '(a, i0, 5a)' )&
   'prtf_ifc is ',dtset%prt_ifc, '. The only allowed values',ch10, &
   'are 0 (no output) or 1 (AI2PS format)',ch10,  &
   'Action: correct prt_ifc in your input file.'
   ABI_ERROR(message)
 end if
! check that ifcout is set
 if (dtset%prt_ifc /= 0 .and. dtset%ifcout == 0) then
   dtset%ifcout = -1  ! this forces output of all IFC
 end if

!Default is no output of the DDB to file
 dtset%prtddb = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'prtddb',tread, 'INT')
 if(tread == 1) dtset%prtddb = intarr(1)
 if(dtset%prtddb < 0 .or. dtset%prtddb > 1) then
   write(message, '(a, i0, 5a)' )&
   'prtf_ddb is ',dtset%prtddb, '. The only allowed values',ch10, &
   'are 0 (no output) or 1 (print DDB and DDB.nc files)',ch10,  &
   'Action: correct prtddb in your input file.'
   ABI_ERROR(message)
 end if
! check that ifcflag is set
 if (dtset%prtddb /= 0 .and. dtset%ifcflag == 0) then
   dtset%ifcflag = 1  ! this forces the use of IFC
 end if

 dtset%prtmbm = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'prtmbm',tread, 'INT')
 if(tread == 1) dtset%prtmbm = intarr(1)
!FIXME: should check whether value of prtmbm is valid

!Default is no output of the nesting factor
 dtset%prtnest = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'prtnest',tread, 'INT')
 if(tread == 1) dtset%prtnest = intarr(1)
 if(dtset%prtnest < 0 .or. dtset%prtnest > 2) then
   write(message, '(a, i0, 5a)' )&
   'prtnest is ',dtset%prtnest, ' The only allowed values',ch10, &
   'are 0 (no nesting), 1 (XY format) or 2 (XY+Xcrysden format)',ch10, &
   'Action: correct prtnest in your input file.'
   ABI_ERROR(message)
 end if

 dtset%prtphbands = 1
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'prtphbands',tread, 'INT')
 if (tread == 1) dtset%prtphbands = intarr(1)
 if (all(dtset%prtphbands /= [0, 1, 2])) then
   write(message, '(a, i0, a)' )&
    'prtphbands is ',dtset%prtphbands, ', but the only allowed values are [0, 1, 2].'
   ABI_ERROR(message)
 end if

 dtset%prtsrlr = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'prtsrlr',tread, 'INT')
 if(tread == 1) dtset%prtsrlr = intarr(1)
 if(dtset%prtsrlr < 0 .or. dtset%prtsrlr > 1)then
   write(message, '(a, i0, 5a)' )&
   'prtsrlr is ',dtset%prtsrlr, ', but the only allowed values',ch10, &
   'are 0 or 1.',ch10, 'Action: correct prtsrlr in your input file.'
   ABI_ERROR(message)
 end if

 dtset%prtvol = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'prtvol',tread, 'INT')
 if(tread == 1) dtset%prtvol = intarr(1)

!Q

 dtset%q2shft(:)=zero
 call intagm(dprarr, intarr, jdtset, marr, 3, string(1:lenstr), 'q2shft',tread, 'DPR')
 if(tread == 1) dtset%q2shft(:)=dprarr(1:3)
!FIXME: need a test on valid entries for q2shft

 dtset%qgrid_type = 1  ! default is uniform nqpt(:) grid
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'qgrid_type',tread, 'INT')
 if(tread == 1) dtset%qgrid_type = intarr(1)
 if(dtset%qgrid_type < 1 .or. dtset%qgrid_type > 2) then
   write(message, '(a, i0, 5a)' )&
   'qgrid_type is ',dtset%qgrid_type, ' The only allowed values',ch10, &
   'are 1 (uniform grid from nqpt) or 2 (listed in ep_nqpt, ep_qptlist)',ch10, &
   'Action: correct qgrid_type in your input file.'
   ABI_ERROR(message)
 end if

 dtset%qrefine = 1  ! default is no refinement
 call intagm(dprarr, intarr, jdtset, marr, 3, string(1:lenstr), 'qrefine',tread, 'INT')
 if(tread == 1) dtset%qrefine = intarr(1:3)
 do ii = 1, 3
   if(dtset%qrefine(ii) < 1) then
     write(message, '(a, 3i0, a, a, a, a, a)' )&
     'qrefine is',dtset%qrefine, ' The only allowed values',ch10, &
     'are integers >= 1 giving the refinement of the ngqpt grid',ch10, &
     'Action: correct qrefine in your input file.'
     ABI_ERROR(message)
   end if
 end do

 dtset%quadquad = 1
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'quadquad',tread, 'INT')
 if(tread == 1) dtset%quadquad = intarr(1)
 if(dtset%quadquad < -1 .or. dtset%quadquad > 1)then
   write(message, '(a, i0, 5a)' )&
   'quadquad is ',dtset%quadquad, ', but the only allowed values',ch10, &
   'are 0 or 1 .',ch10, 'Action: correct quadquad in your input file.'
   ABI_ERROR(message)
 end if

!R

 dtset%ramansr = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'ramansr',tread, 'INT')
 if(tread == 1) dtset%ramansr = intarr(1)
 if(dtset%ramansr < 0 .or. dtset%ramansr > 2)then
   write(message, '(a, i0, 5a)' )&
   'ramansr is ',dtset%ramansr, ', but the only allowed values',ch10, &
   'are 0, 1 or 2.',ch10, 'Action: correct ramansr in your input file.'
   ABI_ERROR(message)
 end if

 dtset%relaxat = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'relaxat',tread, 'INT')
 if(tread == 1) dtset%relaxat = intarr(1)
 if(dtset%relaxat < 0 .or. dtset%relaxat > 1)then
   write(message, '(a, i0, 5a)' )&
   'relaxat is ',dtset%relaxat, ', but the only allowed values',ch10, &
   'are 0 or 1.',ch10, 'Action: correct relaxat in your input file.'
   ABI_ERROR(message)
 end if

 dtset%relaxstr = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'relaxstr',tread, 'INT')
 if(tread == 1) dtset%relaxstr = intarr(1)
 if(dtset%relaxstr < 0 .or. dtset%relaxstr > 1)then
   write(message, '(a, i0, 5a)' )&
   'relaxstr is ',dtset%relaxstr, 'but the only allowed values',ch10, &
   'are 0 or 1.',ch10, 'Action: correct relaxstr in your input file.'
   ABI_ERROR(message)
 end if

 dtset%rfmeth = 1
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'rfmeth',tread, 'INT')
 if(tread == 1) dtset%rfmeth = intarr(1)
 if(dtset%rfmeth < 1 .or. dtset%rfmeth > 2)then
   write(message, '(a, i0, 5a)' )&
   'rfmeth is ',dtset%rfmeth, ', but the only allowed values',ch10, &
   'are 1 or 2.',ch10, 'Action: correct rfmeth in your input file.'
   ABI_ERROR(message)
 end if

 dtset%rifcsph = zero
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'rifcsph',tread, 'DPR')
 if(tread == 1) dtset%rifcsph = dprarr(1)
! if(dtset%rifcsph < -tol12)then
!   write(message, '(a, f10.3, 3a)' )&
!&   'rifcsph is ',dtset%rifcsph, ', which is lower than zero.',ch10, &
!&   'Action: correct rifcsph in your input file.'
!   ABI_ERROR(message)
! end if

!S

 dtset%selectz = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'selectz',tread, 'INT')
 if(tread == 1) dtset%selectz = intarr(1)
 if(dtset%selectz < 0 .or. dtset%selectz > 2)then
   write(message, '(a, i0, 5a)' )&
   'selectz is ',dtset%selectz, ', but the only allowed values',ch10, &
   'are 0, 1 or 2 .',ch10, 'Action: correct selectz in your input file.'
   ABI_ERROR(message)
 end if

 dtset%symdynmat = 1
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'symdynmat',tread, 'INT')
 if(tread == 1) dtset%symdynmat = intarr(1)
 if(dtset%symdynmat /= 0 .and. dtset%symdynmat /= 1)then
   write(message, '(a, i0, 5a)' )&
   'symdynmat is ',dtset%symdynmat, '. The only allowed values',ch10, &
   'are 0, or 1.',ch10, 'Action: correct symdynmat in your input file.'
   ABI_ERROR(message)
 end if

!T

 dtset%targetpol(:) = 0._dp
 call intagm(dprarr, intarr, jdtset, marr, 3, string(1:lenstr), 'targetpol',tread, 'DPR')
 if(tread == 1) dtset%targetpol(1:3) = dprarr(1:3)

!Default is use gaussian integration
 dtset%telphint = 1
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'telphint',tread, 'INT')
 if(tread == 1) dtset%telphint = intarr(1)
 if(dtset%telphint < 0 .or. dtset%telphint > 3) then
   write(message, '(a, i0, 6a)' )&
   'telphint is ',dtset%telphint, '. The only allowed values',ch10, &
   'are 0 (tetrahedron) or 1 (gaussian) or ','2 (set of bands occupied ep_b_min, ep_b_max) or 3 (Fermi Dirac).',ch10, &
   'Action: correct telphint in your input file.'
   ABI_ERROR(message)
 end if

 dtset%temperinc = 100.0_dp
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'temperinc',tread, 'DPR')
 if(tread == 1) dtset%temperinc = dprarr(1)
 if(dtset%temperinc < zero)then
   write(message, '(a, f10.3, 3a)' )&
   'temperinc is ',dtset%temperinc, ', which is lower than 0 .',ch10, &
   'Action: correct temperinc in your input file.'
   ABI_ERROR(message)
 end if

 dtset%tempermin = 100.0_dp
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'tempermin',tread, 'DPR')
 if(tread == 1) dtset%tempermin = dprarr(1)
 if(dtset%tempermin < -tol12)then
   write(message, '(a, f10.3, 3a)' )&
   'tempermin is ',dtset%tempermin, ', which is lower than 0 .',ch10, &
   'Action: correct tempermin in your input file.'
   ABI_ERROR(message)
 end if

 dtset%thermal_supercell(:,:)=0
 marr = 9
 ABI_FREE(intarr)
 ABI_FREE(dprarr)
 ABI_MALLOC(intarr, (marr))
 ABI_MALLOC(dprarr, (marr))
 call intagm(dprarr, intarr, jdtset, marr, 9, string(1:lenstr), 'thermal_supercell',tread, 'INT')
 if(tread == 1) dtset%thermal_supercell(1:3, 1:3)=reshape(intarr(1:9), (/3, 3/))
 call mati3det(dtset%thermal_supercell, idet)
 if(sum(abs(dtset%thermal_supercell))>0 .and. idet == 0) then
   write(message, '(a, 9I6, 5a)' )&
   'thermal_supercell is ',dtset%thermal_supercell, ', but the matrix must be non singular',ch10, &
   'with a non zero determinant.',ch10, 'Action: correct thermal_supercell in your input file.'
   ABI_ERROR(message)
 end if

 dtset%thmflag = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'thmflag',tread, 'INT')
 if(tread == 1) dtset%thmflag = intarr(1)
 if(dtset%thmflag < 0 .or. dtset%thmflag > 8)then
   write(message, '(a, i0, 5a)' )&
   'thmflag is ',dtset%thmflag, ', but the only allowed values',ch10, &
   'are between 0 to 8 (included).',ch10, 'Action: correct thmflag in your input file.'
   ABI_ERROR(message)
 end if

 dtset%thmtol = 0.25_dp
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'thmtol',tread, 'DPR')
 if(tread == 1) dtset%thmtol = dprarr(1)
 if(dtset%thmtol < zero)then
   write(message, '(a, es14.4, 3a)' )&
   'thmtol is ',dtset%thmtol, ', which is lower than 0 .',ch10, &
   'Action: correct thmtol in your input file.'
   ABI_ERROR(message)
 end if

 dtset%ep_prt_yambo = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'ep_prt_yambo',tread, 'INT')
 if(tread == 1) dtset%ep_prt_yambo = intarr(1)
 if(dtset%ep_prt_yambo < 0 .or. dtset%ep_prt_yambo > 1) then
   write(message, '(a, i0, 5a)' )&
   'ep_prt_yambo is ',dtset%ep_prt_yambo, ', but the only allowed values',ch10, &
   'are 0 or 1.',ch10, 'Action: correct ep_prt_yambo in your input file.'
   ABI_ERROR(message)
 end if

!default means _do_ symmetrize the ep coupling matrices over qpoints
 dtset%symgkq = 1
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'symgkq',tread, 'INT')
 if(tread == 1) dtset%symgkq = intarr(1)
 if(dtset%symgkq < 0 .or. dtset%symgkq > 1) then
   write(message, '(a, i0, 5a)' )&
   'symgkq is ',dtset%symgkq, ', but the only allowed values',ch10, &
   'are 0 or 1.',ch10, 'Action: correct symgkq in your input file.'
   ABI_ERROR(message)
 else if (dtset%symgkq == 0) then
   ABI_WARNING('You have turned off el-ph matrix symmetrization over q. Use at own risk')
 end if

!U

 dtset%use_k_fine = 0
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'use_k_fine',tread, 'INT')
 if(tread == 1) dtset%use_k_fine = intarr(1)
 if(dtset%use_k_fine /= 1 .and. dtset%use_k_fine /= 0) then
   write(message, '(a, i0, 5a)' )&
   'use_k_fine is ',dtset%use_k_fine, ', but the only allowed values',ch10, &
   'are 1 or 0.',ch10, 'Action: correct use_k_fine in your input file.'
   ABI_ERROR(message)
 end if

 if(dtset%use_k_fine == 1) then
   if (sum(dtset%kptrlatt) == 0 .or. sum(dtset%kptrlatt_fine) == 0 ) then
     ABI_ERROR('If a finer k-grid is used, you must specify both kptrlatt and kptrlatt_fine')
   end if
 end if


!V
 dtset%vs_qrad_tolkms(:) = zero
 call intagm(dprarr, intarr, jdtset, marr, 2, string(1:lenstr), 'vs_qrad_tolkms',tread, 'DPR')
 if (tread == 1) then
    dtset%vs_qrad_tolkms(:) = dprarr(1:2)
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

 ABI_MALLOC(dtset%atifc, (dtset%natifc))
 ABI_MALLOC(dtset%atifcflg, (natom))
 dtset%atifc(:) = 0
 if(dtset%natifc >= 1)then
   ! default to 1 for first natifc atoms
   dtset%atifc(1:dtset%natifc)=1

   if(dtset%natifc > marr)then
     marr = dtset%natifc
     ABI_FREE(intarr)
     ABI_FREE(dprarr)
     ABI_MALLOC(intarr, (marr))
     ABI_MALLOC(dprarr, (marr))
   end if
   call intagm(dprarr, intarr, jdtset, marr, dtset%natifc, string(1:lenstr), 'atifc',tread, 'INT')
   if(tread == 1) dtset%atifc = intarr(1:dtset%natifc)
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

 ABI_MALLOC(dtset%iatfix, (natom))
 dtset%iatfix(:) = 0
 if ((dtset%relaxat == 1).and.(dtset%natfix > 0)) then
   if(natom > marr)then
     marr = natom
     ABI_FREE(intarr)
     ABI_FREE(dprarr)
     ABI_MALLOC(intarr, (marr))
     ABI_MALLOC(dprarr, (marr))
   end if
   call intagm(dprarr, intarr, jdtset, marr, dtset%natfix, string(1:lenstr), 'iatfix',tread, 'INT')
   if(tread == 1) dtset%iatfix(1:dtset%natfix) = intarr(1:dtset%natfix)
 end if
!FIXME: need a test on values of iatfix: are they just 1 or 0?

 if ((dtset%relaxstr == 1).and.(dtset%nstrfix > 0)) then
   dtset%istrfix(:) = 0
   call intagm(dprarr, intarr, jdtset, marr, dtset%nstrfix, string(1:lenstr), 'istrfix',tread, 'INT')
   if(tread == 1) dtset%istrfix(1:dtset%nstrfix) = intarr(1:dtset%nstrfix)
 end if
!FIXME: need a test on values of istrfix

 if (dtset%natprj_bs > 0) then
   ABI_MALLOC(dtset%iatprj_bs, (dtset%natprj_bs))
   if(dtset%natprj_bs > marr)then
     marr = dtset%natprj_bs
     ABI_FREE(intarr)
     ABI_FREE(dprarr)
     ABI_MALLOC(intarr, (marr))
     ABI_MALLOC(dprarr, (marr))
   end if
   dtset%iatprj_bs(:)=0
   call intagm(dprarr, intarr, jdtset, marr, dtset%natprj_bs, string(1:lenstr), 'iatprj_bs',tread, 'INT')
   if(tread == 1) then
     dtset%iatprj_bs(1:dtset%natprj_bs)=intarr(1:dtset%natprj_bs)
   else
     write(message, '(3a)')&
     'natprj_bs is non zero but iatprj_bs is absent ',ch10, &
     'Action: specify iatprj_bs in your input file.'
     ABI_ERROR(message)
   end if
 end if

!J

!K

!L


 if (dtset%lwfflag .eq. 1 ) then
    ABI_MALLOC(dtset%lwf_anchor_iband, (dtset%lwf_nwann))

    if(dtset%lwf_nwann > marr)then
       marr = dtset%lwf_nwann
       ABI_FREE(intarr)
       ABI_FREE(dprarr)
       ABI_MALLOC(intarr, (marr))
       ABI_MALLOC(dprarr, (marr))
    end if
    dtset%lwf_anchor_iband(:)=0
    call intagm(dprarr, intarr, jdtset, marr, dtset%lwf_nwann, string(1:lenstr), 'lwf_anchor_iband',tread, 'INT')
    if(tread == 1) then
       dtset%lwf_anchor_iband(1:dtset%lwf_nwann)=intarr(1:dtset%lwf_nwann)
    !else
    !   write(message, '(3a)')&
    !        'lwfflag > 0 and lwf_anchor_proj = 1 but lwf_anchor_iband is absent ',ch10, &
    !        'Action: specify lwf_anchor_iband in your input file.'
    !   ABI_ERROR(message)
    end if

 else if (dtset%lwfflag .eq. 2 ) then

    ABI_MALLOC(dtset%lwf_projector, (dtset%lwf_nwann))

    if(dtset%lwf_nwann > marr)then
       marr = dtset%lwf_nwann
       ABI_FREE(intarr)
       ABI_FREE(dprarr)
       ABI_MALLOC(intarr, (marr))
       ABI_MALLOC(dprarr, (marr))
    end if
    dtset%lwf_projector(:)=0
    call intagm(dprarr, intarr, jdtset, marr, dtset%lwf_nwann, string(1:lenstr), 'lwf_projector',tread, 'INT')
    if(tread == 1) then
       dtset%lwf_projector(1:dtset%lwf_nwann)=intarr(1:dtset%lwf_nwann)
    else
       write(message, '(3a)')&
            'lwfflag = 2 and lwf_anchor_proj = 1 but lwf_projector is absent ',ch10, &
            'Action: specify lwf_projector in your input file.'
       ABI_ERROR(message)
    end if

 end if

!M

!N

!O

!P

!Q

 if (dtset%nqshft /= 0)then
   if(3*dtset%nqshft > marr)then
     marr = 3*dtset%nqshft
     ABI_FREE(intarr)
     ABI_FREE(dprarr)
     ABI_MALLOC(intarr, (marr))
     ABI_MALLOC(dprarr, (marr))
   end if
   dtset%q1shft(:,:)=zero
   call intagm(dprarr, intarr, jdtset, marr, 3*dtset%nqshft, string(1:lenstr), 'q1shft',tread, 'DPR')
   if(tread == 1) dtset%q1shft(1:3, 1:dtset%nqshft)=&
&   reshape(dprarr(1:3*dtset%nqshft), (/3, dtset%nqshft/))
 end if

 ABI_MALLOC(dtset%qph1l, (3, dtset%nph1l))
 ABI_MALLOC(dtset%qnrml1, (dtset%nph1l))
 if (dtset%nph1l /= 0)then
   if(4*dtset%nph1l > marr)then
     marr = 4*dtset%nph1l
     ABI_FREE(intarr)
     ABI_FREE(dprarr)
     ABI_MALLOC(intarr, (marr))
     ABI_MALLOC(dprarr, (marr))
   end if
   dtset%qph1l(:,:)=zero
   dtset%qnrml1(:)=zero
   call intagm(dprarr, intarr, jdtset, marr, 4*dtset%nph1l, string(1:lenstr), 'qph1l',tread, 'DPR')
   if(tread == 1)then
     do iph1 = 1, dtset%nph1l
       do ii = 1, 3
         dtset%qph1l(ii, iph1)=dprarr(ii+(iph1-1)*4)
       end do
       dtset%qnrml1(iph1)=dprarr(4+(iph1-1)*4)
       if(abs(dtset%qnrml1(iph1))<DDB_QTOL)then
         write(message, '(5a)' )&
         'The first list of wavevectors ','should not have non-analytical data.',ch10, &
         'Action: correct the first list',' of wavevectors in the input file.'
         ABI_ERROR(message)
       end if
     end do
   end if
 end if

 ABI_MALLOC(dtset%qph2l, (3, dtset%nph2l))
 ABI_MALLOC(dtset%qnrml2, (dtset%nph2l))
 if (dtset%nph2l /= 0)then
   if(4*dtset%nph2l > marr)then
     marr = 4*dtset%nph2l
     ABI_FREE(intarr)
     ABI_FREE(dprarr)
     ABI_MALLOC(intarr, (marr))
     ABI_MALLOC(dprarr, (marr))
   end if
   dtset%qph2l(:,:)=zero
   dtset%qnrml2(:)=zero
   call intagm(dprarr, intarr, jdtset, marr, 4*dtset%nph2l, string(1:lenstr), 'qph2l',tread, 'DPR')
   if(tread == 1)then
     do iph2 = 1, dtset%nph2l
       do ii = 1, 3
         dtset%qph2l(ii, iph2)=dprarr(ii+(iph2-1)*4)
       end do
       dtset%qnrml2(iph2)=dprarr(4+(iph2-1)*4)
       if(abs(dtset%qnrml2(iph2))>DDB_QTOL)then
         write(message, '(5a)' )&
         'The second list of wavevectors',' should have only non-analytical data.',ch10, &
         'Action: correct the second list','of wavevectors in the input file.'
         ABI_ERROR(message)
       end if
     end do
   end if
 end if

 if (dtset%nqpath > 0) then
   ABI_MALLOC(dtset%qpath, (3, dtset%nqpath))
   if(3*dtset%nqpath > marr)then
     marr = 3*dtset%nqpath
     ABI_FREE(intarr)
     ABI_FREE(dprarr)
     ABI_MALLOC(intarr, (marr))
     ABI_MALLOC(dprarr, (marr))
   end if
   dtset%qpath(:,:)=zero
   call intagm(dprarr, intarr, jdtset, marr, 3*dtset%nqpath, string(1:lenstr), 'qpath',tread, 'DPR')
   if(tread == 1) then
     dtset%qpath(1:3, 1:dtset%nqpath)= reshape(dprarr(1:3*dtset%nqpath), (/3, dtset%nqpath/))
   else
     write(message, '(3a)')&
     'nqpath is non zero but qpath is absent ',ch10, &
     'Action: specify qpath in your input file.'
     ABI_ERROR(message)
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
!Finished reading in variables-deallocate
!=======================================================================

 ABI_FREE(dprarr)
 ABI_FREE(intarr)

!=======================================================================
!Check consistency of input variables:
!=======================================================================

 if (dtset%frmin > dtset%frmax) then
   write(message, '(3a)' )&
   'frmax should be higher than frmin',ch10, &
   'Action: change frmax and/or frmin  in your input file.'
   ABI_ERROR(message)
 end if

 if (dtset%nqpath == 0 .and. dtset%elphflag == 1) then
   write(message, '(4a)' )&
   'elphflag is 1 but no nqpath has been specified','for phonon linewidths',ch10, &
   'Action: specify nqpath and qpath(3, nqpath) in your input file.'
   ABI_ERROR(message)
 end if

 if(dtset%telphint /= 2 .and. (dtset%ep_b_min /= 0 .or. dtset%ep_b_max /= 0)) then
   write(message, '(a, i0, 3a)' )&
   'telphint is ',dtset%telphint, ', but ep_b_min or ep_b_max',ch10, &
   'are set /= 1. They will not be used'
   call wrtout(std_out, message, 'COLL')
   ABI_WARNING(message)

 else if(dtset%telphint == 2 .and. (dtset%ep_b_min == 0 .or. dtset%ep_b_max == 0)) then
   write(message, '(a, i0, 4a)' )&
   'telphint is ',dtset%telphint, ', but ep_b_min or ep_b_max',ch10, &
   'are not both set. ',ch10, &
   'Action: set ep_b_min and ep_b_max in your input file.',ch10
   ABI_ERROR(message)
 end if

 if(dtset%thmflag < 3) then
   if ((dtset%telphint == 0 .or. dtset%prtnest == 1 .or. &
        dtset%prtnest == 2 .or. dtset%prtfsurf == 1) .and. sum(dtset%kptrlatt) == 0 ) then
     write (message, '(3a)') &
     'if tetrahedron integration is used, ',&
     'or the output of the nesting function/Fermi surface is required, ',&
     'you must specify the kptrlatt'
     ABI_ERROR(message)
   end if
 end if

 if(dtset%prtdos /= 0 .and. dtset%ifcflag /= 1) then
   write(message, '(3a)' )&
   'ifcflag must be 1 when the calculation of the phonon DOS is required ',ch10, &
   'Action: correct ifcflag in your input file.'
   ABI_ERROR(message)
 end if

 if(dtset%prtsrlr /= 0 .and. dtset%ifcflag /= 1) then
   write(message, '(3a)' )&
   'ifcflag must be 1 for the SR/LR decomposition of the phonon frequencies',ch10, &
   'Action: correct ifcflag in your input file.'
   ABI_ERROR(message)
 end if

 if (dtset%gruns_nddbs /= 0 .and. dtset%ifcflag /= 1) then
   ABI_ERROR("ifcflag must be 1 for Gruneisen calculation")
 end if

 if (dtset%vs_qrad_tolkms(1) /= zero .and. dtset%ifcflag /= 1) then
   ABI_ERROR("ifcflag must be 1 to calculate speed of sound")
 end if

 if(dtset%prtdos /= 0 .and. sum(abs(dtset%ng2qpt(:))) < 3 ) then
   write(message, '(3a)' )&
   'ng2qpt must be specified when the calculation of the phonon DOS is required ',ch10, &
   'Action: correct ng2qpt in your input file.'
   ABI_ERROR(message)
 end if

 if (dtset%ifltransport /= 0 .and. dtset%ep_keepbands /= 1) then
   write(message, '(3a)' )&
   'Band dependency of electron phonon matrix elements must be kept for transport ',ch10, &
   'Action: set ep_keepbands to 1 in your input file.'
   ABI_ERROR(message)
 end if

 if (dtset%ifltransport > 1 .and. sum(abs(dtset%kptrlatt)) == 0) then
   write(message, '(3a)' )&
   'For inelastic transport or electron lifetime calculations you must specify kprtlatt ',ch10, &
   'Action: copy kptrlatt from your abinit GS file to your anaddb input file.'
   ABI_ERROR(message)
 end if

!FIXME: add check that if freeze_displ /= 0 then you need to be doing ifc and phonon interpolation

 if (dtset%ifcflag > 0 .and. sum(abs(dtset%ngqpt)) == 0) then
   write(message, '(3a)' )&
   'if you want interatomic force constant output, anaddb needs ngqpt input variable ',ch10, &
   'Action: set ngqpt in your input file.'
   ABI_ERROR(message)
 end if


 if (dtset%ifcflag /= 1 .and. dtset%lwfflag > 0) then
   write(message, '(3a)' )&
       'if you want to construct the lattice wannier functions, IFC must be computed.', ch10, &
       'Action: set ifcflag to 1.'
   ABI_ERROR(message)
 end if



!check that q-grid refinement is a divisor of ngqpt in each direction
 if(any(dtset%qrefine(1:3) > 1) .and. &
    any(abs(dmod(dble(dtset%ngqpt(1:3))/dble(dtset%qrefine(1:3)), one)) > tol10) ) then
   write(message, '(a, 3i10, a, a, a, 3i8, a, a)' )&
   'qrefine is',dtset%qrefine(1:3), ' The only allowed values',ch10, &
   'are integers which are divisors of the ngqpt grid', dtset%ngqpt(1:3), ch10, &
   'Action: correct qrefine in your input file.'
   ABI_ERROR(message)
 end if

!check that fermie and nelect are not both specified
 if(abs(dtset%elph_fermie) > tol10 .and. abs(dtset%ep_extrael) > tol10) then
   write(message, '(a, E10.2, a, E10.2, a, a, a)' )&
    'elph_fermie (',dtset%elph_fermie, ') and ep_extrael (',dtset%ep_extrael, '), may not both be non 0',ch10, &
    'Action: remove one of the two in your input file.'
   ABI_ERROR(message)
 end if

 ! Check for possible typos.
 call anaddb_chkvars(string)

end subroutine invars9
!!***

!----------------------------------------------------------------------

!!****f*m_anaddb_dataset/outvars_anaddb
!!
!! NAME
!! outvars_anaddb
!!
!! FUNCTION
!! Open input file for the anaddb code, then
!! echoes the input information.
!!
!! INPUTS
!! dtset= (derived datatype) contains all the input variables
!! nunit = unit number for input or output
!!
!! OUTPUT
!!  (only writing)
!!
!! NOTES
!! Should be executed by one processor only.
!!
!! SOURCE

subroutine outvars_anaddb(dtset, nunit)

!Arguments-------------------------------
!scalars
 class(anaddb_dataset_type), intent(inout):: dtset
 integer, intent(in):: nunit

!Local variables-------------------------
!scalars
 integer:: ii, iph1, iph2, iqpt, iqshft

!*********************************************************************

!Write the heading
 write(nunit, '(a, 80a, a)') ch10, ('=',ii = 1, 80), ch10
 write(nunit, '(2a)' )' -outvars_anaddb: echo values of input variables ----------------------',ch10

!The flags
 if (dtset%dieflag /= 0 .or. dtset%ifcflag /= 0 .or. &
     dtset%flexoflag /= 0 .or. &
     dtset%nlflag /= 0 .or. dtset%thmflag /= 0 .or. &
     dtset%elaflag /= 0 .or. dtset%elphflag /= 0 .or. &
     dtset%polflag /= 0 .or. dtset%instrflag /= 0 .or. &
     dtset%piezoflag /= 0) then
   write(nunit, '(a)')' Flags :'
   if(dtset%dieflag /= 0)write(nunit, '(3x, a9, 3i10)')'  dieflag',dtset%dieflag
   if(dtset%flexoflag /= 0)write(nunit, '(3x, a9, 3i10)')'flexoflag',dtset%flexoflag
   if(dtset%ifcflag /= 0)write(nunit, '(3x, a9, 3i10)')'  ifcflag',dtset%ifcflag
   if(dtset%nlflag /= 0)write(nunit, '(3x, a9, 3i10)')'   nlflag',dtset%nlflag
   if(dtset%thmflag /= 0)write(nunit, '(3x, a9, 3i10)')'  thmflag',dtset%thmflag
   if(dtset%elaflag /= 0)write(nunit, '(3x, a9, 3i10)')'  elaflag',dtset%elaflag
   if(dtset%elphflag /= 0)write(nunit, '(3x, a9, 3i10)')' elphflag',dtset%elphflag
   if(dtset%polflag /= 0)write(nunit, '(3x, a9, 3i10)')'  polflag',dtset%polflag
   if(dtset%instrflag /= 0)write(nunit, '(3x, a9, 3i10)')'instrflag',dtset%instrflag
   if(dtset%piezoflag /= 0)write(nunit, '(3x, a9, 3i10)')'piezoflag',dtset%piezoflag
   if(dtset%lwfflag /= 0)write(nunit, '(3x, a9, 3i10)')'lwfflag',dtset%lwfflag
 end if

!Write the general information
 if (dtset%rfmeth /= 1 .or. &
     dtset%enunit /= 0 .or. &
     dtset%eivec /= 0 .or. &
     dtset%asr /= 0 .or. &
     dtset%chneut /= 0 .or. &
     dtset%selectz /= 0 .or. dtset%symdynmat /= 1) then
   write(nunit, '(a)')' Miscellaneous information :'
   if(dtset%rfmeth /= 1)write(nunit, '(3x, a9, 3i10)')'   rfmeth',dtset%rfmeth
   if(dtset%enunit /= 0)write(nunit, '(3x, a9, 3i10)')'   enunit',dtset%enunit
   if(dtset%eivec /= 0) write(nunit, '(3x, a9, 3i10)')'    eivec',dtset%eivec
   if(dtset%asr /= 0)   write(nunit, '(3x, a9, 3i10)')'      asr',dtset%asr
   if(dtset%chneut /= 1)write(nunit, '(3x, a9, 3i10)')'   chneut',dtset%chneut
   if(dtset%selectz /= 0)write(nunit, '(3x, a9, 3i10)')'  selectz',dtset%selectz
   if(dtset%symdynmat /= 1)write(nunit, '(3x, a9, 3i10)')'symdynmat',dtset%symdynmat
 end if
 if(dtset%prtvol /= 0) write(nunit, '(3x, a9, i10)')'   prtvol',dtset%prtvol

!Frequency information
 if(dtset%dieflag == 1)then
   write(nunit, '(a)')' Frequency information :'
   write(nunit, '(3x, a9, 3i10)')'    nfreq',dtset%nfreq
   write(nunit, '(3x, a9, 7x, 3es16.8)')'    frmin',dtset%frmin
   write(nunit, '(3x, a9, 7x, 3es16.8)')'    frmax',dtset%frmax
 end if

!For interatomic force constant information
 if(dtset%ifcflag /= 0)then
   write(nunit, '(a)')' Interatomic Force Constants Inputs :'
   write(nunit, '(3x, a9, 3i10)')'   dipdip',dtset%dipdip
   write(nunit, '(3x, a9, 3i10)')'   dipquad',dtset%dipquad
   write(nunit, '(3x, a9, 3i10)')'   quadquad',dtset%quadquad
   if(dtset%nsphere /= 0)write(nunit, '(3x, a9, 3i10)')'  nsphere',dtset%nsphere
   if(abs(dtset%rifcsph)>tol10)write(nunit, '(3x, a9, E16.6)')'  nsphere',dtset%rifcsph
   write(nunit, '(3x, a9, 3i10)')'   ifcana',dtset%ifcana
   write(nunit, '(3x, a9, 3i10)')'   ifcout',dtset%ifcout
   if(dtset%natifc >= 1)then
     write(nunit, '(3x, a9, 3i10)')'   natifc',dtset%natifc
     write(nunit, '(3x, a9, 8i10)')'    atifc',(dtset%atifc(ii), ii = 1, dtset%natifc)
   end if
   write(nunit, '(a)')' Description of grid 1 :'
   write(nunit, '(3x, a9, 3i10)')'     brav',dtset%brav
   write(nunit, '(3x, a9, 3i10)')'    ngqpt',dtset%ngqpt(1:3)
   write(nunit, '(3x, a9, 3i10)')'   nqshft',dtset%nqshft
   if (dtset%nqshft /= 0)then
     write(nunit, '(3x, a9)')'   q1shft'
     do iqshft = 1, dtset%nqshft
       write(nunit, '(19x, 4es16.8)') (dtset%q1shft(ii, iqshft), ii = 1, 3)
     end do
   end if
   if (any(dtset%qrefine(:) > 1)) then
     write(nunit, '(3x, a9, 3i10)')'  qrefine', dtset%qrefine
   end if
   ! Speed of sound
   if (dtset%vs_qrad_tolkms(1) > zero) then
      write(nunit, '(a, 2es16.8)')"vs_qrad_tolkms", (dtset%vs_qrad_tolkms(:))
   end if
 end if

!Phonon density of states with gaussian method
 if(dtset%prtdos /= 0)then
   write(nunit, '(a)')' Phonon DOS information :'
   write(nunit, '(3x, a9, es16.8)')'dosdeltae',dtset%dosdeltae
   write(nunit, '(3x, a9, es16.8)')' dossmear',dtset%dossmear
 end if

!Thermal information
 if(dtset%thmflag /= 0)then
   write(nunit, '(a)')' Thermal information :'
   write(nunit, '(3x, a9, 3i10)')'    nchan',dtset%nchan
   write(nunit, '(3x, a9, 3i10)')'   nwchan',dtset%nwchan
   write(nunit, '(3x, a9, 7x, 3es16.8)')'   dostol',dtset%dostol
   write(nunit, '(3x, a9, 7x, 3es16.8)')'   thmtol',dtset%thmtol
   write(nunit, '(3x, a9, 3i10)')'  ntemper',dtset%ntemper
   write(nunit, '(3x, a9, 7x, 3es16.8)')'temperinc',dtset%temperinc
   write(nunit, '(3x, a9, 7x, 3es16.8)')'tempermin',dtset%tempermin
 endif

!Grid 2 description
 if(dtset%thmflag /= 0 .or. dtset%prtdos /= 0)then
   write(nunit, '(a)')' Description of grid 2 (Fourier interp. or BZ sampling):'
   write(nunit, '(3x, a9, 3i10)')'   ng2qpt',dtset%ng2qpt(1:3)
   write(nunit, '(3x, a9, 3i10)')'   ngrids',dtset%ngrids
   write(nunit, '(3x, a9, 7x, 3es16.8)')'   q2shft',dtset%q2shft(1:3)
 end if

!Non-linear response information
 if (dtset%nlflag /= 0) then
   write(nunit, '(a)')' Non-linear response information :'
   write(nunit, '(3x, a9, i10)') '   alphon',dtset%alphon
   write(nunit, '(3x, a9, 3i10)')'   prtmbm',dtset%prtmbm
   write(nunit, '(3x, a9, 3i10)')'  ramansr',dtset%ramansr
 end if

!Structural relaxation at fixed polarization
 if (dtset%polflag /= 0) then
   write(nunit, '(a)')' Relaxation at fixed polarization :'
   if (dtset%relaxat == 1) then
     write(nunit, '(3x, a9, i10)') '  relaxat',dtset%relaxat
   end if
   if (dtset%relaxstr == 1) then
     write(nunit, '(a12, i10)') ' relaxstr',dtset%relaxstr
   end if
 end if

!Elphon information
 if (dtset%elphflag /= 0) then
   write(nunit, '(a)')' Elphon calculation will be carried out'
   write(nunit, '(a12, E16.6)') 'elphsmear', dtset%elphsmear
   write(nunit, '(a12, E16.6)') 'a2fsmear', dtset%a2fsmear
   write(nunit, '(a12, E16.6)') 'mustar', dtset%mustar
   write(nunit, '(a12, i10)') 'nqpath', dtset%nqpath
   write(nunit, '(a12)') 'qpath'
   do iqpt = 1, dtset%nqpath
     write(nunit, '(12x, 3(E16.6, 1x))') dtset%qpath(:,iqpt)
   end do
   write(nunit, '(a12, i10)') 'telphint', dtset%telphint
   if (dtset%telphint == 0) then
     write(nunit, '(a)') ' Tetrahedron integration for elphon'
   else if (dtset%telphint == 1) then
     write(nunit, '(a)') ' Smeared weight integration for elphon'
   else if (dtset%telphint == 2) then
     write(nunit, '(a)') ' Band filtered integration for elphon'
   end if
   if (abs(dtset%elph_fermie) > tol10) then
     write(nunit, '(a12, E16.6)')  'elph_fermie', dtset%elph_fermie
   end if
   if (dtset%ep_extrael /= 0) then
     if (abs(dtset%ep_extrael) > 1.0d2) then
        write(nunit, '(a, E20.12)')' Doping set by the user is (negative for el doping) :',dtset%ep_extrael
     else
       write(nunit, '(a, E16.6)')  'Elphon: extra electrons per unit cell = ', dtset%ep_extrael
     end if
   end if
   if (dtset%ep_nspline /= 20) then
     write(nunit, '(a, I8)')  'Elphon: scale factor for spline interpolation in RTA = ', dtset%ep_nspline
   end if
   if (dtset%band_gap < 10.0d0) then
     write(nunit, '(a, E16.6)')  'Elphon: set band gap to (in eV) = ', dtset%band_gap
   end if

   if (sum(abs(dtset%kptrlatt)) > 0) then
     write(nunit, '(a12, 3(3(i3, 1x), 2x))' ) 'kptrlatt',reshape( dtset%kptrlatt(:,:), (/9/) )
   end if

   if (sum(abs(dtset%kptrlatt_fine)) > 0) then
     write(nunit, '(a12, 3(3(i3, 1x), 2x))' ) 'kptrlatt_fine ',reshape( dtset%kptrlatt_fine(:,:), (/9/) )
   end if

   if (dtset%ep_keepbands == 1) then
     write(nunit, '(a)') ' Will keep band dependency in gkk in memory.'
     write(nunit, '(a)') ' WARNING: the memory requirements will be multiplied by nbands**2 !!!'
   end if

   if (dtset%ep_scalprod == 1) then
     write(nunit, '(a)') ' scalar product will be performed when assembling the gamma matrices.'
     write(nunit, '(a)') ' WARNING: with this option you can not distinguish which '
     write(nunit, '(a)') '    linewidth comes from which phonon mode !!!'
   end if

   if (dtset%prtbltztrp == 1) write(nunit, '(a)') ' Will output input files for BoltzTraP'
   if (dtset%prtfsurf == 1) write(nunit, '(a)') ' Will output fermi surface in XCrysDen format'
   if (dtset%prt_ifc == 1) write(nunit, '(a)') ' Will output real space IFC in AI2PS and TDEP format'
   if (dtset%prtnest == 1) write(nunit, '(a)') ' Will output nesting factor'

   if (dtset%ifltransport == 1) then
     write(nunit, '(a)') ' Will perform transport calculation in elphon to get'
     write(nunit, '(a, a)') ' resistivity and thermal conductivity as a function of T',ch10
     write(nunit, '(a, es16.6, a)' ) ' Minimum temperature for transport outputs: ', dtset%tempermin, ' K'
     write(nunit, '(a, es16.6, a)' ) ' Maximum temperature for transport outputs: ', &
       dtset%tempermin+dtset%temperinc*dtset%ntemper, ' K'
     write(nunit, '(a, i6)' ) ' Number of temperature points for transport outputs: ', dtset%ntemper
     write(nunit, '(a)' )
   end if

   if (dtset%gkqwrite == 1) then
     write(nunit, '(a, a)' ) 'Gkk matrix elements on input grid of ',&
     'qpoints will be written to disk. File gkqfile must be absent.'
   end if
   if (dtset%gkk_rptwrite == 1) then
     write(nunit, '(a, a)' ) 'Gkk matrix elements in real space ',&
     'will be written to disk. File gkk_rpt_file must be absent.'
   end if
   if (dtset%gkk2write == 1) then
     write(nunit, '(a, a)' ) 'Full grid gkk matrix elements ',&
     'will be written to disk. File gkk2file must be absent.'
   end if
 end if

 if (dtset%gruns_nddbs /= 0) then
   write(nunit, '(a)' ) "Will compute Gruneisen parameters with finite difference method. DDB files:"
   do ii = 1, dtset%gruns_nddbs
     write(nunit, "(2a)")"    ",trim(dtset%gruns_ddbs(ii))
   end do
 end if

! lattice wannier function Information
 if (dtset%lwfflag > 0) then
    write(nunit, '(a)')' Lattice Wannier function information:'
    write(nunit, '(a20, i10)')    '            lwfflag', dtset%lwfflag
    write(nunit, '(a20, i10)')    '          lwf_nwann', dtset%lwf_nwann
    write(nunit, '(a20, i10)')    '    lwf_anchor_proj', dtset%lwf_anchor_proj
    write(nunit, '(a20, 3i10)')   '          lwf_ngqpt',(dtset%lwf_ngqpt(ii), ii = 1, 3)
    write(nunit, '(a20, i10)')    '    lwf_disentangle', dtset%lwf_disentangle
    write(nunit, '(a20, E16.6)')  '             lwf_mu', dtset%lwf_mu
    write(nunit, '(a20, E16.6)')  '          lwf_sigma', dtset%lwf_sigma
    write(nunit, '(a20, 3E16.6)') '     lwf_anchor_qpt',(dtset%lwf_anchor_qpt(ii), ii = 1, 3)
    if (abs(dtset%lwf_anchor_proj) > 0) then
        write(nunit, '(a20)',advance="no")      '   lwf_anchor_iband'
        do ii = 1, dtset%lwf_nwann
          write(nunit, '(3x, I5)', advance="no") dtset%lwf_anchor_iband(ii)
        end do
        write(nunit, '(a)') ' '
    end if

    if (dtset%lwfflag .eq. 2) then
       write(nunit, '(a20)',advance="no")      '    lwf_projector'
       do ii = 1, dtset%lwf_nwann
          write(nunit, '(3x, I5)', advance="no") dtset%lwf_projector(ii)
       end do
       write(nunit, '(a)') ' '
    end if
 end if


!List of vector 1  (reduced coordinates)
 if(dtset%nph1l /= 0)then
   write(nunit, '(a)')' First list of wavevector (reduced coord.) :'
   write(nunit, '(3x, a9, 3i10)')'    nph1l',dtset%nph1l
   write(nunit, '(3x, a9)')'    qph1l'
   do iph1 = 1, dtset%nph1l
     write(nunit, '(19x, 3es16.8, 2x, es11.3)') &
       (dtset%qph1l(ii, iph1), ii = 1, 3), dtset%qnrml1(iph1)
   end do
 end if

!List of vector 2  (cartesian coordinates)
 if(dtset%nph2l /= 0)then
   write(nunit, '(a)')' Second list of wavevector (cart. coord.) :'
   write(nunit, '(3x, a9, 3i10)')'    nph2l',dtset%nph2l
   write(nunit, '(3x, a9)')'    qph2l'
   do iph2 = 1, dtset%nph2l
     write(nunit, '(19x, 3es16.8, 2x, es11.3)') (dtset%qph2l(ii, iph2), ii = 1, 3), dtset%qnrml2(iph2)
   end do
 end if

!phonon frozen in supercell
 if (abs(dtset%freeze_displ) > tol10) then
   write(nunit, '(a)') 'Phonon displacements will be output, frozen into supercells'
   write(nunit, '(a, E20.10)') ' Chosen amplitude of frozen displacements = ', dtset%freeze_displ
 end if

!atom projected bs files
 if (abs(dtset%natprj_bs) > 0) then
   write(nunit, '(a)') 'Phonon band structure files, with atomic projections, will be output '
   write(nunit, '(a)') ' Chosen atoms for projection = '
   write(nunit, '(10I6)') dtset%iatprj_bs
 end if

 write(nunit, '(a, 80a, a)') ch10, ('=',ii = 1, 80), ch10

end subroutine outvars_anaddb
!!***

!----------------------------------------------------------------------

!!****f*m_anaddb_dataset/anaddb_dtset_init
!!
!! NAME
!! anaddb_dtset_init
!!
!! FUNCTION
!! Initialize the code ppddb9: write heading and make the first i/os
!!
!! INPUTS
!!  input_path: String with input file path. Empty string activates files file in legacy mode.
!!
!! OUTPUT
!! character(len = fnlen) filnam(7)=character strings giving file names
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
!! SOURCE

subroutine anaddb_dtset_init(dtset, input_path)

!Arguments-------------------------------
 class(anaddb_dataset_type), intent(inout):: dtset
 character(len=*), intent(in):: input_path

!Local variables-------------------------
!scalars
 integer:: lenstr, marr, jdtset, tread, i1, ierr
 character(len=strlen):: string, raw_string, fname, dirpath
 character(len=fnlen):: unused
 character(len = 500):: msg
!arrays
 integer, allocatable:: intarr(:)
 real(dp), allocatable:: dprarr(:)

! *********************************************************************

 dtset%filename_input = input_path
 dtset%filename_output = "run.abo"
 dtset%prefix_outdata = trim(" ")
 dtset%filename_ddb = trim(" ")
 dtset%filename_ddk = trim(" ")
 dtset%prefix_eph = trim(" ")
 dtset%filename_gkk = trim(" ")
 dtset%filename_eigr2d = trim(" ")


 if (len_trim(input_path) == 0) then
   !write(msg, "(3a)") "Please run Anaddb as: anaddb input.in",ch10,&
   !                   "Note that using a files file is no longer supported in Abinit10."
   !ABI_ERROR(msg)
   ! Legacy Files file mode.
   write(std_out, "(2a)")" DeprecationWarning: ",ch10
   write(std_out, "(a)") "     The files file has been deprecated in Abinit9 and will be removed in Abinit10."
   write(std_out, "(2a)")"     Use the syntax `anaddb t01.abi` where t01.abi is an anaddb input with ddb_filepath.",ch10
   write(std_out, "(3a)")'            ddb_filepath = "out_DDB"',ch10, ch10

   write(std_out, *)' Give name for formatted input file: '
   read(std_in, '(a)' ) dtset%filename_input
   write(std_out, '(a, a)' )'-   ',trim(dtset%filename_input)
   write(std_out, *)' Give name for formatted output file: '
   read(std_in, '(a)' ) dtset%filename_output
   write(std_out, '(a, a)' )'-   ',trim(dtset%filename_output)
   write(std_out, *)' Give name for input derivative database: '
   read(std_in, '(a)' ) dtset%filename_ddb
   write(std_out, '(a, a)' )'-   ',trim(dtset%filename_ddb)
   write(std_out, *)' Give name for output molecular dynamics: '
   read(std_in, '(a)' ) unused
   write(std_out, '(a, a)' )'-   ',trim(unused)
   ! GA: This message is confusing, because filnam(5) is also for EIG2D files.
   write(std_out, *)' Give name for input elphon matrix elements (GKK file): '
   read(std_in, '(a)' ) dtset%filename_gkk
   write(std_out, '(a, a)' )'-   ',trim(dtset%filename_gkk)
   write(std_out, *)' Give root name for elphon output files: '
   read(std_in, '(a)' ) dtset%prefix_eph
   write(std_out, '(a, a)' )'-   ',trim(dtset%prefix_eph)
   write(std_out, *)' Give name for file containing ddk filenames for elphon/transport: '
   read(std_in, '(a)' ) dtset%filename_ddk
   write(std_out, '(a, a)' )'-   ',trim(dtset%filename_ddk)
   dtset%prefix_outdata = trim(" ")

 end if

 ! Read input
 string = repeat(" ", strlen)
 raw_string = repeat(" ", strlen)
 call instrng(dtset%filename_input, lenstr, 1, strlen, string, raw_string)
 ! To make case-insensitive, map characters to upper case.
 call inupper(string(1:lenstr))

 marr = 3
 ABI_MALLOC(intarr, (marr))
 ABI_MALLOC(dprarr, (marr))
 jdtset = 0

 ! Allow user to override default values
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), "output_file", tread, 'KEY', key_value=dtset%filename_output)
 write(std_out, '(2a)')'- Name for formatted output file: ', trim(dtset%filename_output)

 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'ddb_filepath', tread, 'KEY', key_value=dtset%filename_ddb)
 !Check that we didnt use a files file
 if (len_trim(input_path) > 0) then
   ABI_CHECK(tread == 1, 'ddb_filepath variable must be specified in the input file')
 end if
 write(std_out, "(2a)")'- Input derivative database: ', trim(dtset%filename_ddb)

 ! Nobody knows the scope of this line in the files file.
 !call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), 'md_output', tread, 'KEY', key_value=filnam(4))

 ! GA: This variable name is confusing, because filnam(5) is also for EIG2D files.
 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), "gkk_filepath", tread, 'KEY', key_value=dtset%filename_gkk)
 if (tread == 1) write(std_out, "(2a)")'- Name for input elphon matrix elements (GKK file): ', trim(dtset%filename_gkk)
 ! GA: This is to keep old behavior. To be cleaned or removed.
 dtset%filename_eigr2d = dtset%filename_gkk

 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), "eph_prefix", tread, 'KEY', key_value=dtset%prefix_eph)
 if (tread == 1) write(std_out, "(2a)")"- Root name for elphon output files: ", trim(dtset%prefix_eph)

 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), "ddk_filepath", tread, 'KEY', key_value=dtset%filename_ddk)
 if (tread == 1) write(std_out, "(2a)")"- File containing ddk filenames for elphon/transport: ", trim(dtset%filename_ddk)

 call intagm(dprarr, intarr, jdtset, marr, 1, string(1:lenstr), "outdata_prefix", tread, 'KEY', key_value=dtset%prefix_outdata)
 if (tread == 1) then
   write(std_out, "(2a)")'- Root name for output files: ', trim(dtset%prefix_outdata)
 end if

 ABI_FREE(intarr)
 ABI_FREE(dprarr)

 ! Compute OUTPUT_PREFIX as in abinit.
 ! I do not change the "files" file to avoid backward compatibility issue
 if (len_trim(dtset%prefix_outdata) == 0) then
   fname = basename(trim(dtset%filename_output))
   i1 = index(fname, ".",back=.true.)
   if ( i1 > 1 ) then
     dtset%prefix_outdata = fname(:i1-1)
   end if
   write(std_out, "(2a)")'- Root name for output files set to: ', trim(dtset%prefix_outdata)
 endif

 i1 = index(dtset%prefix_outdata, "/", back=.True.)
 if (i1 > 0) then
   dirpath = dtset%prefix_outdata(1:i1-1)
   call clib_mkdir_if_needed(dirpath, ierr)
   ABI_CHECK(ierr == 0, sjoin("Error", itoa(ierr), "while trying to create directory", dirpath))
 end if

end subroutine anaddb_dtset_init
!!***

!!****f*m_anaddb_dataset/anaddb_dtset_bcast_files
!! NAME
!! anaddb_dtset_bcast_files
!!
!! FUNCTION
!!  Broadcast filenames
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine anaddb_dtset_bcast_files(dtset, comm)

 class(anaddb_dataset_type), intent(inout):: dtset
 integer, intent(in) :: comm

 integer, parameter :: master = 0
 integer :: ierr

 call xmpi_bcast(dtset%filename_input, master, comm, ierr)
 call xmpi_bcast(dtset%filename_output, master, comm, ierr)
 call xmpi_bcast(dtset%prefix_outdata, master, comm, ierr)
 call xmpi_bcast(dtset%filename_ddb, master, comm, ierr)
 call xmpi_bcast(dtset%filename_ddk, master, comm, ierr)
 call xmpi_bcast(dtset%prefix_eph, master, comm, ierr)
 call xmpi_bcast(dtset%filename_gkk, master, comm, ierr)
 call xmpi_bcast(dtset%filename_eigr2d, master, comm, ierr)

end subroutine anaddb_dtset_bcast_files
!!***

!!****f*m_anaddb_dataset/anaddb_chkvars
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
!! SOURCE

subroutine anaddb_chkvars(string)

!Arguments------------------------------------
!scalars
 character(len=*), intent(in):: string

!Local variables-------------------------------
!scalars
 integer, parameter:: protocol0 = 0
 character(len = 100):: list_logicals, list_strings, list_vars_img
 character(len = 10000):: list_vars

!************************************************************************

!Here, list all admitted variable names (max 10 per line, to fix the ideas)
!Note: Do not use "double quotation mark" for the string since it triggers a bug in docchk.py (abirules script)
!<ANADDB_VARS>
!A
 list_vars=                 ' alphon asr a2fsmear atifc'
!B
 list_vars = trim(list_vars)//' brav band_gap'
!C
 list_vars = trim(list_vars)//' chneut'
!D
 list_vars = trim(list_vars)//' dieflag dipdip dipquad dossum dosdeltae dossmear dostol dos_maxmode'
!E
 list_vars = trim(list_vars)//' ep_scalprod eivec elaflag elphflag enunit'
 list_vars = trim(list_vars)//' ep_b_min ep_b_max ep_int_gkk ep_keepbands ep_nqpt ep_nspline ep_prt_yambo'
 list_vars = trim(list_vars)//' elphsmear elph_fermie ep_extrael ep_qptlist'
!F
 list_vars = trim(list_vars)//' flexoflag freeze_displ frmax frmin'
!G
 list_vars = trim(list_vars)//' gkk2write gkk_rptwrite gkqwrite gruns_nddbs'
!H
!I
 list_vars = trim(list_vars)//' ifcana ifcflag ifcout ifltransport instrflag istrfix iatfix iatprj_bs'
!J
!K
 list_vars = trim(list_vars)//' kptrlatt kptrlatt_fine'
!L

 list_vars = trim(list_vars)//' lwf_anchor_iband lwf_anchor_proj lwf_anchor_qpt'
 list_vars = trim(list_vars)//' lwf_disentangle lwf_mu lwf_ngqpt lwf_nwann lwf_projector lwf_sigma'
 list_vars = trim(list_vars)//' lwfflag'
!M
 list_vars = trim(list_vars)//' mustar'
!N
 list_vars = trim(list_vars)//' natfix natifc natom natprj_bs nchan ndivsm nfreq ngrids nlflag nph1l nph2l'
 list_vars = trim(list_vars)//' nqpath nqshft nsphere nstrfix ntemper nwchan ngqpt ng2qpt'
!O
 list_vars = trim(list_vars)//' outboltztrap'
!P
 list_vars = trim(list_vars)//' piezoflag polflag prtddb prtdos prt_ifc prtmbm prtfsurf'
 list_vars = trim(list_vars)//' prtnest prtphbands prtsrlr prtvol prtbltztrp'
!Q
 list_vars = trim(list_vars)//' qrefine qgrid_type q1shft q2shft qnrml1 qnrml2 qpath qph1l qph2l quadquad'
!R
 list_vars = trim(list_vars)//' ramansr relaxat relaxstr rfmeth rifcsph'
!S
 list_vars = trim(list_vars)//' selectz symdynmat symgkq'
!T
 list_vars = trim(list_vars)//' targetpol telphint thmflag temperinc tempermin thermal_supercell thmtol'
!U
 list_vars = trim(list_vars)//' use_k_fine'
!V
 list_vars = trim(list_vars)//' vs_qrad_tolkms'
!W
!X
!Y
!Z

!
 list_vars_img=' '

!Logical input variables
 list_logicals=' '

!String input variables
 list_strings=' gruns_ddbs ddb_filepath output_file outdata_prefix gkk_filepath eph_prefix ddk_filepath' ! md_output
!</ANADDB_VARS>

!Extra token, also admitted:
!<ANADDB_UNITS>
 list_vars = trim(list_vars)//' au Angstr Angstrom Angstroms Bohr Bohrs eV meV Ha'
 list_vars = trim(list_vars)//' Hartree Hartrees K nm Ry Rydberg Rydbergs S Sec Second T Tesla'
!</ANADDB_UNITS>

!<ANADDB_OPERATORS>
 list_vars = trim(list_vars)//' sqrt '
!</ANADDB_OPERATORS>

!Transform to upper case
 call inupper(list_vars)
 call inupper(list_vars_img)
 call inupper(list_logicals)
 call inupper(list_strings)

 call chkvars_in_string(protocol0, list_vars, list_vars_img, list_logicals, list_strings, string)

end subroutine anaddb_chkvars
!!***


!!****f*m_anaddb_dataset/anaddb_dtset_read_input
!! NAME
!! anaddb_dtset_read_input
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine anaddb_dtset_read_input(dtset, comm)

 class(anaddb_dataset_type), intent(inout):: dtset
 integer, intent(in):: comm

!Local variables-------------------------------
 integer, parameter:: master = 0
 integer :: lenstr
 integer:: ierr
 integer :: my_rank
 logical :: iam_master
 character(len = strlen):: string, raw_string
 type(ddb_hdr_type) :: ddb_hdr

 my_rank = xmpi_comm_rank(comm)
 iam_master = (my_rank == master)

 ! Must read natom from the DDB before being able to allocate some arrays needed for invars9
 call ddb_hdr%open_read(dtset%filename_ddb, comm=comm, dimonly = 1)

 dtset%natom = ddb_hdr%natom
 dtset%ntypat = ddb_hdr%ntypat
 dtset%mpert = ddb_hdr%mpert
 dtset%msize = ddb_hdr%msize
 dtset%usepaw = ddb_hdr%usepaw

 ! Read the input file, and store the information in a long string of characters
 ! strlen from defs_basis module
 if (iam_master) then
   call instrng(dtset%filename_input, lenstr, 1, strlen, string, raw_string)
   ! To make case-insensitive, map characters to upper case.
   call inupper(string(1:lenstr))
 end if

 call xmpi_bcast(string, master, comm, ierr)
 call xmpi_bcast(raw_string, master, comm, ierr)
 call xmpi_bcast(lenstr, master, comm, ierr)

 dtset%lenstr = lenstr
 dtset%input_string = string

 ! Save input string in global variable so that we can access it in ntck_open_create
 ABI_MALLOC_TYPE_SCALAR(character(len=len_trim(raw_string)), INPUT_STRING)
 INPUT_STRING = trim(raw_string)

 ! Read the inputs
 call invars9(dtset, dtset%lenstr, dtset%natom, dtset%input_string)

 ! Set some inputs depending on what the ddb contains
 !if (.not. ddb_hdr%has_d3E_lw) then
 !  ! The default value is 1.
 !  ! Here we set the flags to zero if Q*is not available.
 !  ! GA: I think this check is correct, but originally it was done
 !  ! by checking for a non-zero value of 
 !  ! iblock_quadrupoles = ddb_lw%get_quadrupoles(ddb_hdr%ddb_version, lwsym, BLKTYP_d3E_lw, qdrp_cart)
 !  dtset%dipquad = 0
 !  dtset%quadquad = 0
 !end if

 call ddb_hdr%free()

 ! GA: This used to be after the dry run exit.
 ! Check the value and transform the meaning of atifc (1 and 0 only)
 call chkin9(dtset%atifcflg,dtset%atifc,dtset%natifc,dtset%natom)

end subroutine anaddb_dtset_read_input
!!***

end module m_anaddb_dataset
!!***

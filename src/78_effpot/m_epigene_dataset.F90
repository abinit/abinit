!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_epigene_dataset
!! NAME
!!  m_epigene_dataset
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2015 ABINIT group (AM)
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

module m_epigene_dataset
    
 use defs_basis
 use m_profiling_abi
 use m_errors

 use m_ddb,  only : DDB_QTOL

 implicit none

 private

 public :: epigene_dataset_type
 public :: epigene_dtset_free
 public :: outvars_epigene
 public :: invars10
!!***

!----------------------------------------------------------------------

!!****t* m_epigene_dataset/epigene_dataset_type
!! NAME
!! epigene_dataset_type
!!
!! FUNCTION
!! The epigene_dataset_type structured datatype
!! gather all the input variables for the epigene code.
!!
!! SOURCE

 type epigene_dataset_type

! Variables should be declared on separated lines in order to reduce the occurence of bzr conflicts.
! Since all these input variables are described in the epigene_help.html
! file, they are not described in length here ...
! Integer

  integer :: asr
  integer :: brav
  integer :: chneut
  integer :: dipdip
  integer :: eivec
  integer :: elphflag
  integer :: enunit
  integer :: ifcana
  integer :: ifcflag
  integer :: ifcout
  integer :: dynamics
  integer :: natifc
  integer :: natom
  integer :: ntime
  integer :: nph1l
  integer :: nph2l
  integer :: nqshft
  integer :: nsphere
  integer :: prt_effpot
  integer :: prt_phfrq
  integer :: prt_ifc
  integer :: prt_3rd  ! Print the 3rd order in xml file
  integer :: prtsrlr  ! print the short-range/long-range decomposition of phonon freq.
  integer :: qrefine
  integer :: rfmeth
  integer :: symdynmat

  integer :: n_cell(3)
  integer :: ngqpt(9)             ! ngqpt(9) instead of ngqpt(3) is needed in wght9.f
  integer :: ng2qpt(3)
  integer :: kptrlatt(3,3)
  integer :: kptrlatt_fine(3,3)

! Real(dp)
  real(dp) :: temperature
  real(dp) :: delta_df
  real(dp) :: energy_reference
  real(dp) :: rifcsph

  real(dp) :: acell(3)
  real(dp) :: strain(6)
  real(dp) :: rprim(3,3)
  real(dp) :: q1shft(3,4)

! Integer arrays
  integer, allocatable :: atifc(:)
   ! atifc(natom)

! Real arrays 
  real(dp), allocatable :: qnrml1(:) 
  ! qnrml1(nph1l)

  real(dp), allocatable :: qnrml2(:)
  ! qnrml1(nph1l)

  real(dp), allocatable :: qph1l(:,:) 
  ! qph1l(3,nph1l)

  real(dp), allocatable :: qph2l(:,:)
  ! qph2l(3,nph2l)

 end type epigene_dataset_type
!!***

contains 
!!***

!!****f* m_epigene_dataset/epigene_dtset_free
!!
!! NAME
!!   epigene_dtset_free
!!
!! FUNCTION
!!   deallocate remaining arrays in the epigene_dtset datastructure
!!
!! INPUTS
!!  epigene_dtset = epigene datastructure
!!
!! PARENTS
!!      epigene
!!
!! CHILDREN
!!
!! NOTES
!!
!! SOURCE

subroutine epigene_dtset_free(epigene_dtset)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'epigene_dtset_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(epigene_dataset_type), intent(inout) :: epigene_dtset

! *************************************************************************
 
 if (allocated(epigene_dtset%atifc))  then
   ABI_DEALLOCATE(epigene_dtset%atifc)
 end if
 if (allocated(epigene_dtset%qnrml1))  then
   ABI_DEALLOCATE(epigene_dtset%qnrml1)
 end if
 if (allocated(epigene_dtset%qnrml2))  then
   ABI_DEALLOCATE(epigene_dtset%qnrml2)
 end if
 if (allocated(epigene_dtset%qph1l))  then
   ABI_DEALLOCATE(epigene_dtset%qph1l)
 end if
 if (allocated(epigene_dtset%qph2l))  then
   ABI_DEALLOCATE(epigene_dtset%qph2l)
 end if

end subroutine epigene_dtset_free
!!***

!----------------------------------------------------------------------

!!****f* m_epigene_dataset/invars10
!!
!! NAME
!! invars9
!!
!! FUNCTION
!! Open input file for the epigene code, then reads or echoes the input information.
!!
!! INPUTS
!! lenstr=actual length of string
!! natom=number of atoms, needed for atifc
!! string*(*)=string of characters containing all input variables and data
!!
!! OUTPUT
!! epigene_dtset= (derived datatype) contains all the input variables
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
!!      epigene
!!
!! CHILDREN
!!
!! SOURCE

subroutine invars10(epigene_dtset,lenstr,natom,string)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'invars10'
 use interfaces_42_parser
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: lenstr,natom
 character(len=*),intent(in) :: string
 type(epigene_dataset_type),intent(inout) :: epigene_dtset 

!Local variables -------------------------
!Dummy arguments for subroutine 'intagm' to parse input file
!Set routine version number here:
!scalars
 integer,parameter :: vrsddb=100401
 integer :: iatifc,ii,iph1,iph2,jdtset,marr,tread
 character(len=500) :: message
!arrays
 integer,allocatable :: intarr(:)
 real(dp),allocatable :: dprarr(:),work(:)

!*********************************************************************
 marr=3
 ABI_ALLOCATE(intarr,(marr))
 ABI_ALLOCATE(dprarr,(marr))

 jdtset=1

!copy natom to epigene_dtset
 epigene_dtset%natom=natom

!=====================================================================
!start reading in dimensions and non-dependent variables
!=====================================================================

!A
 epigene_dtset%asr=2
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'asr',tread,'INT')
 if(tread==1) epigene_dtset%asr=intarr(1)
 if(epigene_dtset%asr<-2.or.epigene_dtset%asr>5)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'asr is',epigene_dtset%asr,', but the only allowed values',ch10,&
&   'are 0, 1, 2, 3, 4, 5, -1 or -2 .',ch10,&
&   'Action: correct asr in your input file.'
!  Note : negative values are allowed when the acoustic sum rule
!  is to be applied after the analysis of IFCs
!  3,4 are for rotational invariance (under development)
!  5 is for hermitian imposition of the ASR
   MSG_ERROR(message)
 end if

!B
 epigene_dtset%brav=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'brav',tread,'INT')
 if(tread==1) epigene_dtset%brav=intarr(1)
 if(epigene_dtset%brav/=1)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'brav is',epigene_dtset%brav,', but the only allowed values',ch10,&
&   'are 1 for epigene (not implemented) .',ch10,&
&   'Action: correct brav in your input file.'
   MSG_ERROR(message)
 end if

!C
 epigene_dtset%chneut=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'chneut',tread,'INT')
 if(tread==1) epigene_dtset%chneut=intarr(1)
 if(epigene_dtset%chneut<0.or.epigene_dtset%chneut>2)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'chneut is',epigene_dtset%chneut,', but the only allowed values',ch10,&
&   'are 0, 1 or 2 .',ch10,&
&   'Action: correct chneut in your input file.'
   MSG_ERROR(message)
 end if

!D
 epigene_dtset%dipdip=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dipdip',tread,'INT')
 if(tread==1) epigene_dtset%dipdip=intarr(1)
 if(epigene_dtset%dipdip>1.or.epigene_dtset%dipdip<0)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'dipdip is',epigene_dtset%dipdip,', but the only allowed values',ch10,&
&   'is 1.',ch10,&
&   'Action: correct dipdip in your input file.'
   MSG_ERROR(message)
 end if

 epigene_dtset%delta_df= 1d-02
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'delta_df',tread,'DPR')
 if(tread==1) epigene_dtset%delta_df=dprarr(1)
 if(epigene_dtset%delta_df<0)then
   write(message, '(a,es10.2,a,a,a,a,a)' )&
&   'delta_df is',epigene_dtset%delta_df,', but the only allowed values',ch10,&
&   'are superior to 0  .',ch10,&
&   'Action: correct delta_df in your input file.'
   MSG_ERROR(message)
 end if

 epigene_dtset%energy_reference= zero
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'energy_reference',tread,'DPR')
 if(tread==1) epigene_dtset%energy_reference=dprarr(1)

 epigene_dtset%enunit=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'enunit',tread,'INT')
 if(tread==1) epigene_dtset%enunit=intarr(1)
 if(epigene_dtset%enunit<0.or.epigene_dtset%enunit>2)then
   write(message, '(a,i0,a,a,a,a,a)' )&
&   'enunit is',epigene_dtset%enunit,', but the only allowed values',ch10,&
&   'are 0, 1 or 2.',ch10,&
&   'Action: correct enunit in your input file.'
   MSG_ERROR(message)
 end if

 epigene_dtset%ifcana=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ifcana',tread,'INT')
 if(tread==1) epigene_dtset%ifcana=intarr(1)
 if(epigene_dtset%ifcana<0.or.epigene_dtset%ifcana>1)then
   write(message, '(a,i0,a,a,a,a,a)' )&
&   'ifcana is',epigene_dtset%ifcana,', but the only allowed values',ch10,&
&   'are 0 or 1.',ch10,&
&   'Action: correct ifcana in your input file.'
   MSG_ERROR(message)
 end if

 epigene_dtset%ifcflag=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ifcflag',tread,'INT')
 if(tread==1) epigene_dtset%ifcflag=intarr(1)
 if(epigene_dtset%ifcflag<0.or.epigene_dtset%ifcflag>1)then
   write(message, '(a,i0,a,a,a,a,a)' )&
&   'ifcflag is',epigene_dtset%ifcflag,', but the only allowed values',ch10,&
&   'are 0 or 1.',ch10,&
&   'Action: correct ifcflag in your input file.'
   MSG_ERROR(message)
 end if

 epigene_dtset%prtsrlr=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtsrlr',tread,'INT')
 if(tread==1) epigene_dtset%prtsrlr=intarr(1)
 if(epigene_dtset%prtsrlr<0.or.epigene_dtset%prtsrlr>1)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'prtsrlr is',epigene_dtset%prtsrlr,', but the only allowed values',ch10,&
&   'are 0 or 1.',ch10,&
&   'Action: correct prtsrlr in your input file.'
   MSG_ERROR(message)
 end if

 epigene_dtset%ifcout=2000000
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ifcout',tread,'INT')
 if(tread==1) epigene_dtset%ifcout=intarr(1)
 if(epigene_dtset%ifcout<0)then
   write(message, '(a,i0,a,a,a)' )&
&   'ifcout is',epigene_dtset%ifcout,', which is lower than 0 .',ch10,&
&   'Action: correct ifcout in your input file.'
   MSG_ERROR(message)
 end if

 epigene_dtset%ntime=200
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ntime',tread,'INT')
 if(tread==1) epigene_dtset%ntime=intarr(1)
 if(epigene_dtset%ntime<0)then
   write(message, '(a,i0,a,a,a)' )&
&   'ntime is',epigene_dtset%ntime,', which is lower than 0 .',ch10,&
&   'Action: correct ntime in your input file.'
   MSG_ERROR(message)
 end if


 epigene_dtset%dynamics=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dynamics',tread,'INT')
 if(tread==1) epigene_dtset%dynamics=intarr(1)
 if(epigene_dtset%dynamics<0.or.epigene_dtset%dynamics>2)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'dynamics is',epigene_dtset%dynamics,', but the only allowed values',ch10,&
&   'are 0 or  1.',ch10,&
&   'Action: correct dynamics in your input file.'
   MSG_ERROR(message)
 end if

!N
 epigene_dtset%natifc=natom
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'natifc',tread,'INT')
 if(tread==1) epigene_dtset%natifc=intarr(1)
 if(epigene_dtset%natifc<0)then
   write(message, '(a,i0,a,a,a)' )&
&   'natifc is',epigene_dtset%natifc,', which is lower than 0 .',ch10,&
&   'Action: correct natifc in your input file.'
   MSG_ERROR(message)
 end if

 if(epigene_dtset%natifc>natom)then
   write(message, '(a,i0,a,a,a,i0,a,a,a)' )&
&   'The number of atom ifc in the input files',epigene_dtset%natifc,',',ch10,&
&   'is larger than the number of atoms',natom,'.',ch10,&
&   'Action: change natifc in the input file.'
   MSG_ERROR(message)
 end if

 epigene_dtset%ng2qpt(:)=0
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'ng2qpt',tread,'INT')
 if(tread==1) epigene_dtset%ng2qpt(:)=intarr(1:3)
 do ii=1,3
   if(epigene_dtset%ng2qpt(ii)<0)then
     write(message, '(a,i0,a,i0,a,a,a,i0,a)' )&
&     'ng2qpt(',ii,') is',epigene_dtset%ng2qpt(ii),', which is lower than 0 .',ch10,&
&     'Action: correct ng2qpt(',ii,') in your input file.'
     MSG_ERROR(message)
   end if
 end do

 epigene_dtset%n_cell(:)= one
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'n_cell',tread,'INT')
 if(tread==1) epigene_dtset%n_cell(1:3)=intarr(1:3)
 do ii=1,3
   if(epigene_dtset%n_cell(ii)<0.or.epigene_dtset%n_cell(ii)>20)then
     write(message, '(a,i0,a,i0,a,a,a,i0,a)' )&
&     'n_cell(',ii,') is ',epigene_dtset%n_cell(ii),', which is lower than 0 of superior than 10.',&
&     ch10,'Action: correct n_cell(',ii,') in your input file.'
     MSG_ERROR(message)
   end if
 end do

 epigene_dtset%ngqpt(:)= one
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'ngqpt',tread,'INT')
 if(tread==1) epigene_dtset%ngqpt(1:3)=intarr(1:3)
 do ii=1,3
   if(epigene_dtset%ngqpt(ii)<0)then
     write(message, '(a,i0,a,i0,a,a,a,i0,a)' )&
&     'ngqpt(',ii,') is',epigene_dtset%ngqpt(ii),', which is lower than 0 .',ch10,&
&     'Action: correct ngqpt(',ii,') in your input file.'
     MSG_ERROR(message)
   end if
 end do

 epigene_dtset%nph1l=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nph1l',tread,'INT')
 if(tread==1) epigene_dtset%nph1l=intarr(1)
 if(epigene_dtset%nph1l<0)then
   write(message, '(a,i0,a,a,a)' )&
&   'nph1l is',epigene_dtset%nph1l,', which is lower than 0 .',ch10,&
&   'Action: correct nph1l in your input file.'
   MSG_ERROR(message)
 end if
 
 epigene_dtset%nph2l=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nph2l',tread,'INT')
 if(tread==1) epigene_dtset%nph2l=intarr(1)
 if(epigene_dtset%nph2l<0)then
   write(message, '(a,i0,a,a,a)' )&
&   'nph2l is',epigene_dtset%nph2l,', which is lower than 0 .',ch10,&
&   'Action: correct nph2l in your input file.'
   MSG_ERROR(message)
 end if
 
 epigene_dtset%nqshft=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nqshft',tread,'INT')
 if(tread==1) epigene_dtset%nqshft=intarr(1)
 if(epigene_dtset%nqshft<0 .or. epigene_dtset%nqshft==3 .or.&
& epigene_dtset%nqshft>=5 )then
   write(message, '(a,i0,a,a,a,a,a)' )&
&   'nqshft is',epigene_dtset%nqshft,', but the only allowed values',ch10,&
&   'are 1, 2 or 4 .',ch10,&
&   'Action: correct nqshft in your input file.'
   MSG_ERROR(message)
 end if

 epigene_dtset%nsphere=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nsphere',tread,'INT')
 if(tread==1) epigene_dtset%nsphere=intarr(1)
 if(epigene_dtset%nsphere<0)then
   write(message, '(a,i0,a,a,a)' )&
&   'nsphere is',epigene_dtset%nsphere,', which is lower than 0',ch10,&
&   'Action: correct nsphere in your input file.'
   MSG_ERROR(message)
 end if

!O

!P
 epigene_dtset%prt_effpot=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prt_effpot',tread,'INT')
 if(tread==1) epigene_dtset%prt_effpot=intarr(1)
 if(epigene_dtset%prt_effpot<-2.or.epigene_dtset%prt_effpot>3)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'prt_effpot is',epigene_dtset%prtsrlr,', but the only allowed values',ch10,&
&   'are 0, 1 or 2.',ch10,&
&   'Action: correct prt_effpot in your input file.'
   MSG_ERROR(message)
 end if

 epigene_dtset%prt_phfrq=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prt_phfrq',tread,'INT')
 if(tread==1) epigene_dtset%prt_phfrq=intarr(1)
 if(epigene_dtset%prt_phfrq<0.or.epigene_dtset%prt_phfrq>2)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'prt_phfrq is',epigene_dtset%prtsrlr,', but the only allowed values',ch10,&
&   'are 0, 1 or 2.',ch10,&
&   'Action: correct prt_phfrq in your input file.'
   MSG_ERROR(message)
 end if

!Default is no output of the real space IFC to file
 epigene_dtset%prt_ifc = 0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prt_ifc',tread,'INT')
 if(tread==1) epigene_dtset%prt_ifc = intarr(1)
 if(epigene_dtset%prt_ifc < 0 .or. epigene_dtset%prt_ifc > 1) then
   write(message, '(a,i0,a,a,a,a,a)' )&
&   'prtf_ifc is',epigene_dtset%prt_ifc,'. The only allowed values',ch10,&
&   'are 0 (no output) or 1 (AI2PS format)',ch10,  &
&   'Action: correct prt_ifc in your input file.'
   MSG_ERROR(message)
 end if

!Default is no output of the 3rd derivative
 epigene_dtset%prt_3rd = 0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prt_3rd',tread,'INT')
 if(tread==1) epigene_dtset%prt_3rd = intarr(1)
 if(epigene_dtset%prt_3rd < 0 .or. epigene_dtset%prt_3rd > 2) then
   write(message, '(a,i0,a,a,a,a,a,a,a)' )&
&   'prtf_3rd is ',epigene_dtset%prt_3rd,'. The only allowed values',ch10,&
&   'are 0 (no computation), 1 (only computation)',ch10,&
&   'or 2 (computation and print in xml file)',ch10,  &
&   'Action: correct prt_3rd in your input file.'
   MSG_ERROR(message)
 end if

!Q
 epigene_dtset%qrefine=1 ! default is no refinement
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'qrefine',tread,'INT')
 if(tread==1) epigene_dtset%qrefine = intarr(1)
 if(epigene_dtset%qrefine < 1) then
   write(message, '(a,i0,a,a,a,a,a)' )&
&   'qrefine is',epigene_dtset%qrefine,' The only allowed values',ch10,&
&   'are integers >= 1 giving the refinement of the ngqpt grid',ch10,&
&   'Action: correct qrefine in your input file.'
   MSG_ERROR(message)
 end if

!R
 epigene_dtset%rfmeth=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rfmeth',tread,'INT')
 if(tread==1) epigene_dtset%rfmeth=intarr(1)
 if(epigene_dtset%rfmeth<1.or.epigene_dtset%rfmeth>2)then
   write(message, '(a,i0,a,a,a,a,a)' )&
&   'rfmeth is',epigene_dtset%rfmeth,', but the only allowed values',ch10,&
&   'are 1 or 2 . ',ch10,&
&   'Action: correct rfmeth in your input file.'
   MSG_ERROR(message)
 end if

 epigene_dtset%rifcsph=zero
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rifcsph',tread,'DPR')
 if(tread==1) epigene_dtset%rifcsph=dprarr(1)
 if(epigene_dtset%rifcsph<-tol12)then
   write(message, '(a,f10.3,a,a,a)' )&
&   'rifcsph is',epigene_dtset%rifcsph,', which is lower than zero.',ch10,&
&   'Action: correct rifcsph in your input file.'
   MSG_ERROR(message)
 end if

!S
 epigene_dtset%symdynmat=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'symdynmat',tread,'INT')
 if(tread==1) epigene_dtset%symdynmat=intarr(1)
 if(epigene_dtset%symdynmat/=0.and.epigene_dtset%symdynmat/=1)then
   write(message, '(a,i0,a,a,a,a,a)' )&
&   'symdynmat is',epigene_dtset%symdynmat,'. The only allowed values',ch10,&
&   'are 0, or 1.',ch10,&
&   'Action: correct symdynmat in your input file.'
   MSG_ERROR(message)
 end if

!T
 epigene_dtset%temperature=325
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'temperature',tread,'DPR')
 if(tread==1) epigene_dtset%temperature=dprarr(1)
 if(epigene_dtset%temperature<=0)then
   write(message, '(a,f10.1,a,a,a,a,a)' )&
&   'Temperature is ',epigene_dtset%temperature,'. The only allowed values',ch10,&
&   'are positives values.',ch10,&
&   'Action: correct Temperature in your input file.'
   MSG_ERROR(message)
 end if

!U

!V

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
 epigene_dtset%acell= one
 if(3>marr)then
   marr=3
   ABI_DEALLOCATE(intarr)
   ABI_DEALLOCATE(dprarr)
   ABI_ALLOCATE(intarr,(marr))
   ABI_ALLOCATE(dprarr,(marr))
 end if
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'acell',tread,'DPR')
 if(tread==1) epigene_dtset%acell(1:3)= dprarr(1:3)
 if(any(epigene_dtset%acell<=tol10))then
    write(message, '(3a)' )&
&       'There is negative on zero value for cell ',ch10,&
&       'Action: change acell in your input file.'
      MSG_ERROR(message) 
 end if

 ABI_ALLOCATE(epigene_dtset%atifc,(natom))
 epigene_dtset%atifc(:)=zero
 if(epigene_dtset%natifc>=1)then
   if(epigene_dtset%natifc>marr)then
     marr=epigene_dtset%natifc
     ABI_DEALLOCATE(intarr)
     ABI_DEALLOCATE(dprarr)
     ABI_ALLOCATE(intarr,(marr))
     ABI_ALLOCATE(dprarr,(marr))
   end if
   call intagm(dprarr,intarr,jdtset,marr,epigene_dtset%natifc,string(1:lenstr),'atifc',tread,'INT')
   if(tread==1) then 
     epigene_dtset%atifc(1:epigene_dtset%natifc)= intarr(1:epigene_dtset%natifc)
   else ! set to the maximum
     do iatifc=1,epigene_dtset%natifc
       epigene_dtset%atifc(iatifc) =  iatifc
     end do
   end if
   ABI_MALLOC(work,(natom))
   work(:)=0

   do iatifc=1,epigene_dtset%natifc
     if(epigene_dtset%atifc(iatifc)<=0.or.epigene_dtset%atifc(iatifc)>natom)then
       write(message, '(a,i0,a,a,a,a,a,i0,a,a,a)' )&
&       'For iatifc=',iatifc,', the number of the atom ifc to be ',ch10,&
&       'analysed is not valid : either negative, ',ch10,&
&       'zero, or larger than natom =',natom,'.',ch10,&
&       'Action: change atifc in your input file.'
       MSG_ERROR(message)
     end if
     work(epigene_dtset%atifc(iatifc))=1
   end do
   epigene_dtset%atifc(1:natom)=work(:)
   ABI_FREE(work)
 end if

!B

!C

!D

!E
 epigene_dtset%eivec=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'eivec',tread,'INT')
 if(tread==1) epigene_dtset%eivec=intarr(1)
 if(epigene_dtset%eivec<0.or.epigene_dtset%eivec>4)then
   write(message, '(a,i0,a,a,a,a,a)' )&
&   'eivec is',epigene_dtset%eivec,', but the only allowed values',ch10,&
&   'are 0, 1, 2, 3 or 4.',ch10,&
&   'Action: correct eivec in your input file.'
   MSG_ERROR(message)
 end if

!F

!G

!H

!I


!J

!K

!L

!M

!N

!O

!P

!Q

 if (epigene_dtset%nqshft/=0)then
   if(3*epigene_dtset%nqshft>marr)then
     marr=3*epigene_dtset%nqshft
     ABI_DEALLOCATE(intarr)
     ABI_DEALLOCATE(dprarr)
     ABI_ALLOCATE(intarr,(marr))
     ABI_ALLOCATE(dprarr,(marr))
   end if
   epigene_dtset%q1shft(:,:)=zero
   call intagm(dprarr,intarr,jdtset,marr,3*epigene_dtset%nqshft, string(1:lenstr),'q1shft',tread,'DPR')
   if(tread==1) epigene_dtset%q1shft(1:3,1:epigene_dtset%nqshft)=&
&   reshape(dprarr(1:3*epigene_dtset%nqshft),(/3,epigene_dtset%nqshft/))
 end if

 ABI_ALLOCATE(epigene_dtset%qph1l,(3,epigene_dtset%nph1l))
 ABI_ALLOCATE(epigene_dtset%qnrml1,(epigene_dtset%nph1l))
 if (epigene_dtset%nph1l/=0)then
   if(4*epigene_dtset%nph1l>marr)then
     marr=4*epigene_dtset%nph1l
     ABI_DEALLOCATE(intarr)
     ABI_DEALLOCATE(dprarr)
     ABI_ALLOCATE(intarr,(marr))
     ABI_ALLOCATE(dprarr,(marr))
   end if
   epigene_dtset%qph1l(:,:)=zero
   epigene_dtset%qnrml1(:)=zero
   call intagm(dprarr,intarr,jdtset,marr,4*epigene_dtset%nph1l,string(1:lenstr),'qph1l',tread,'DPR')
   if(tread==1)then
     do iph1=1,epigene_dtset%nph1l
       do ii=1,3
         epigene_dtset%qph1l(ii,iph1)=dprarr(ii+(iph1-1)*4)
       end do
       epigene_dtset%qnrml1(iph1)=dprarr(4+(iph1-1)*4)
       if(abs(epigene_dtset%qnrml1(iph1))<DDB_QTOL)then
         write(message, '(a,a,a,a,a)' )&
&         'The first list of wavevectors ','should not have non-analytical data.',ch10,&
&         'Action: correct the first list',' of wavevectors in the input file.'
         MSG_ERROR(message)
       end if
     end do
   end if
 end if

 ABI_ALLOCATE(epigene_dtset%qph2l,(3,epigene_dtset%nph2l))
 ABI_ALLOCATE(epigene_dtset%qnrml2,(epigene_dtset%nph2l))
 if (epigene_dtset%nph2l/=0)then
   if(4*epigene_dtset%nph2l>marr)then
     marr=4*epigene_dtset%nph2l
     ABI_DEALLOCATE(intarr)
     ABI_DEALLOCATE(dprarr)
     ABI_ALLOCATE(intarr,(marr))
     ABI_ALLOCATE(dprarr,(marr))
   end if
   epigene_dtset%qph2l(:,:)=zero
   epigene_dtset%qnrml2(:)=zero
   call intagm(dprarr,intarr,jdtset,marr,4*epigene_dtset%nph2l,string(1:lenstr),'qph2l',tread,'DPR')
   if(tread==1)then
     do iph2=1,epigene_dtset%nph2l
       do ii=1,3
         epigene_dtset%qph2l(ii,iph2)=dprarr(ii+(iph2-1)*4)
       end do
       epigene_dtset%qnrml2(iph2)=dprarr(4+(iph2-1)*4)
       if(abs(epigene_dtset%qnrml2(iph2))>DDB_QTOL)then
         write(message, '(a,a,a,a,a)' )&
&         'The second list of wavevectors',' should have only non-analytical data.',ch10,&
&         'Action: correct the second list','of wavevectors in the input file.'
         MSG_ERROR(message)
       end if
     end do
   end if
 end if

!R
 if(9>marr)then
   marr=9
   ABI_DEALLOCATE(intarr)
   ABI_DEALLOCATE(dprarr)
   ABI_ALLOCATE(intarr,(marr))
   ABI_ALLOCATE(dprarr,(marr))
 end if
 epigene_dtset%rprim(:,:)= zero
 call intagm(dprarr,intarr,jdtset,marr,9,string(1:lenstr),'rprim',tread,'DPR')
 if(tread==1) then
   epigene_dtset%rprim(1:3,1:3)= reshape(dprarr(1:9),(/3,3/))
! check new rprimd
   if(all(epigene_dtset%rprim(1,:)==zero).or.&
&    all(epigene_dtset%rprim(2,:)==zero).or.all(epigene_dtset%rprim(3,:)==zero)) then
     write(message, '(3a)' )&
&  ' There is a problem with rprim',ch10,&
&   'Action: correct rprim'
     MSG_BUG(message)
   end if
 end if
!S

 if(6>marr)then
   marr=6
   ABI_DEALLOCATE(intarr)
   ABI_DEALLOCATE(dprarr)
   ABI_ALLOCATE(intarr,(marr))
   ABI_ALLOCATE(dprarr,(marr))
 end if
 epigene_dtset%strain(:)= zero
 call intagm(dprarr,intarr,jdtset,marr,6,string(1:lenstr),'strain',tread,'DPR')
 if(tread==1) epigene_dtset%strain(1:6)= dprarr(1:6)

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

 if(epigene_dtset%prtsrlr/=0 .and. epigene_dtset%ifcflag/=1) then
   write(message, '(3a)' )&
&   'ifcflag must be 1 for the SR/LR decomposition of the phonon frequencies',ch10,&
&   'Action: correct ifcflag in your input file.'
   MSG_ERROR(message)
 end if

!FIXME: add check that if freeze_displ /= 0 then you need to be doing ifc and phonon interpolation

 if (epigene_dtset%ifcflag > 0 .and. sum(abs(epigene_dtset%ngqpt)) == 0) then
   write(message, '(3a)' )&
&   'if you want interatomic force constant output, epigene needs ngqpt input variable ',ch10,&
&   'Action: set ngqpt in your input file.'
   MSG_ERROR(message)
 end if

!check that q-grid refinement is a divisor of ngqpt in each direction
 if(epigene_dtset%qrefine > 1 .and. sum(abs(dmod(epigene_dtset%ngqpt/dble(epigene_dtset%qrefine),one))) > tol10) then
   write(message, '(a,i0,a,a,a,3i8,a,a)' )&
&   'qrefine is',epigene_dtset%qrefine,' The only allowed values',ch10,&
&   'are integers which are divisors of the ngqpt grid', epigene_dtset%ngqpt,ch10,&
&   'Action: correct qrefine in your input file.'
   MSG_ERROR(message)
 end if

! check new rprimd
 if(all(epigene_dtset%acell(:) > one).and.all(epigene_dtset%rprim(:,:)==zero))then
   write(message, '(3a)' )&
&         ' acell is defined but there is no rprim',ch10,&
&         'Action: add rprim input'
   MSG_BUG(message)
 end if


end subroutine invars10
!!***

!----------------------------------------------------------------------

!!****f* m_epigene_dataset/outvars_epigene
!!
!! NAME
!! outvars_epigene
!!
!! FUNCTION
!! Open input file for the epigene code, then
!! echoes the input information.
!!
!! INPUTS
!! epigene_dtset= (derived datatype) contains all the input variables
!! nunit=unit number for input or output
!!
!! OUTPUT
!!  (only writing)
!!
!! NOTES
!! Should be executed by one processor only.
!!
!! PARENTS
!!      epigene
!!
!! CHILDREN
!!
!! SOURCE

subroutine outvars_epigene (epigene_dtset,nunit)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'outvars_epigene'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: nunit
 type(epigene_dataset_type),intent(in) :: epigene_dtset

!Local variables -------------------------
!Set routine version number here:
!scalars
 integer :: ii,iph1,iph2,iqshft

!*********************************************************************

!Write the heading
 write(nunit,'(a,80a,a)') ch10,('=',ii=1,80),ch10
 write(nunit, '(a,a)' )&
& ' -outvars_epigene: echo values of input variables ----------------------',ch10

!The flags
 if(epigene_dtset%ifcflag/=0)then
   write(nunit,'(a)')' Flags :'
   if(epigene_dtset%ifcflag/=0)write(nunit,'(3x,a9,3i10)')'  ifcflag',epigene_dtset%ifcflag
   if(epigene_dtset%prt_effpot/=0)write(nunit,'(3x,a9,3i10)')'prt_effpot',epigene_dtset%prt_effpot
   if(epigene_dtset%prt_phfrq/=0)write(nunit,'(3x,a9,3i10)')'prt_phfrq',epigene_dtset%prt_phfrq
   if(epigene_dtset%prt_3rd/=0)write(nunit,'(3x,a9,3i10)')'  prt_3rd',epigene_dtset%prt_3rd
   if(epigene_dtset%prt_3rd==2)write(nunit,'(3x,a9,3es8.2)')'delta_df',epigene_dtset%delta_df
 end if

 if(epigene_dtset%dynamics/=0)then
   write(nunit,'(a)')' Molecular Dynamics :'
   write(nunit,'(3x,a9,3F10.1)')'     temp',epigene_dtset%temperature
   write(nunit,'(3x,a9,3I10.1)')'    ntime',epigene_dtset%ntime
   write(nunit,'(3x,a9,3i10)')  '    ncell',epigene_dtset%n_cell
 end if

!Write the general information
 if( epigene_dtset%rfmeth/=1 .or. &
& epigene_dtset%enunit/=0 .or. &
& epigene_dtset%eivec/=0 .or. &
& epigene_dtset%asr/=0 .or. &
& epigene_dtset%chneut/=0)then
   write(nunit,'(a)')' Miscellaneous information :'
   if(epigene_dtset%rfmeth/=1)write(nunit,'(3x,a9,3i10)')'   rfmeth',epigene_dtset%rfmeth
   if(epigene_dtset%enunit/=0)write(nunit,'(3x,a9,3i10)')'   enunit',epigene_dtset%enunit
   if(epigene_dtset%eivec/=0) write(nunit,'(3x,a9,3i10)')'    eivec',epigene_dtset%eivec
   if(epigene_dtset%asr/=0)   write(nunit,'(3x,a9,3i10)')'      asr',epigene_dtset%asr
   if(epigene_dtset%chneut/=0)write(nunit,'(3x,a9,3i10)')'   chneut',epigene_dtset%chneut
 end if


!For interatomic force constant information
 if(epigene_dtset%ifcflag/=0)then
   write(nunit,'(a)')' Interatomic Force Constants Inputs :'
   write(nunit,'(3x,a9,3i10)')'   dipdip',epigene_dtset%dipdip
   if(epigene_dtset%nsphere/=0)write(nunit,'(3x,a9,3i10)')'  nsphere',epigene_dtset%nsphere
   if(abs(epigene_dtset%rifcsph)>tol10)write(nunit,'(3x,a9,E16.6)')'  nsphere',epigene_dtset%rifcsph
   write(nunit,'(3x,a9,3i10)')'   ifcana',epigene_dtset%ifcana
   write(nunit,'(3x,a9,3i10)')'   ifcout',epigene_dtset%ifcout
   if(epigene_dtset%natifc>=1)then
     write(nunit,'(3x,a9,3i10)')'   natifc',epigene_dtset%natifc
     write(nunit,'(3x,a12)',advance='no')'    atifc   '
     write(nunit,'(3x,15i4)') (epigene_dtset%atifc(ii)*ii,ii=1,epigene_dtset%natifc)

   end if
   write(nunit,'(a)')' Description of grid 1 :'
   write(nunit,'(3x,a9,3i10)')'     brav',epigene_dtset%brav
   write(nunit,'(3x,a9,3i10)')'    ngqpt',epigene_dtset%ngqpt(1:3)
   write(nunit,'(3x,a9,3i10)')'   nqshft',epigene_dtset%nqshft
   if (epigene_dtset%nqshft/=0)then
     write(nunit,'(3x,a9)')'   q1shft'
     do iqshft=1,epigene_dtset%nqshft
       write(nunit,'(19x,4es16.8)') (epigene_dtset%q1shft(ii,iqshft),ii=1,3)
     end do
   end if
   if (epigene_dtset%qrefine > 1) then
     write(nunit,'(3x,a9,i10)')'  qrefine', epigene_dtset%qrefine
   end if
 end if


!List of vector 1  (reduced coordinates)
 if(epigene_dtset%nph1l/=0)then
   write(nunit,'(a)')' First list of wavevector (reduced coord.) :'
   write(nunit,'(3x,a9,3i10)')'    nph1l',epigene_dtset%nph1l
   write(nunit,'(3x,a9)')'    qph1l'
   do iph1=1,epigene_dtset%nph1l
     write(nunit,'(19x,3es16.8,2x,es11.3)') &
&     (epigene_dtset%qph1l(ii,iph1),ii=1,3),epigene_dtset%qnrml1(iph1)
   end do
 end if

!List of vector 2  (cartesian coordinates)
 if(epigene_dtset%nph2l/=0)then
   write(nunit,'(a)')' Second list of wavevector (cart. coord.) :'
   write(nunit,'(3x,a9,3i10)')'    nph2l',epigene_dtset%nph2l
   write(nunit,'(3x,a9)')'    qph2l'
   do iph2=1,epigene_dtset%nph2l
     write(nunit,'(19x,3es16.8,2x,es11.3)') &
&     (epigene_dtset%qph2l(ii,iph2),ii=1,3),epigene_dtset%qnrml2(iph2)
   end do
 end if

 write(nunit,'(a,80a,a)') ch10,('=',ii=1,80),ch10

end subroutine outvars_epigene
!!***

!----------------------------------------------------------------------

end module m_epigene_dataset
!!***

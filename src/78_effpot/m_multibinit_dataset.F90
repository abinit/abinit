!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_multibinit_dataset
!! NAME
!!  m_multibinit_dataset
!!
!! FUNCTION
!!  module with the type for the input of multibinit
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2017 ABINIT group (AM)
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

module m_multibinit_dataset
    
 use defs_basis
 use m_profiling_abi
 use m_errors

 use m_ddb,  only : DDB_QTOL

 implicit none

 private

 public :: multibinit_dataset_type
 public :: multibinit_dtset_free
 public :: outvars_multibinit
 public :: invars10
!!***

!----------------------------------------------------------------------

!!****t* m_multibinit_dataset/multibinit_dataset_type
!! NAME
!! multibinit_dataset_type
!!
!! FUNCTION
!! The multibinit_dataset_type structured datatype
!! gather all the input variables for the multibinit code.
!!
!! SOURCE

 type multibinit_dataset_type

! Variables should be declared on separated lines in order to reduce the occurence of bzr conflicts.
! Since all these input variables are described in the multibinit_help.html
! file, they are not described in length here ...
! Integer
  integer :: asr
  integer :: brav
  integer :: chneut
  integer :: confinement
  integer :: conf_power_disp
  integer :: conf_power_strain
  integer :: dipdip
  integer :: eivec
  integer :: elphflag
  integer :: enunit
  integer :: fit_coeff
  integer :: fit_option
  integer :: fit_ncycle
  integer :: fit_nfixcoeff
  integer :: ifcana
  integer :: ifcflag
  integer :: ifcout
  integer :: dtion
  integer :: dynamics
  integer :: natifc
  integer :: natom
  integer :: ncoeff
  integer :: ntime
  integer :: nnos
  integer :: nph1l
  integer :: nph2l
  integer :: nqshft
  integer :: nsphere
  integer :: optcell
  integer :: prt_model
  integer :: prt_phfrq
  integer :: prt_ifc
  integer :: strcpling  ! Print the 3rd order in xml file
  integer :: prtsrlr  ! print the short-range/long-range decomposition of phonon freq.
  integer :: rfmeth
  integer :: restarxf
  integer :: symdynmat

  integer :: fit_grid(3)
  integer :: fit_rangePower(2)
  integer :: n_cell(3)
  integer :: ngqpt(9)             ! ngqpt(9) instead of ngqpt(3) is needed in wght9.f
  integer :: ng2qpt(3)
  integer :: kptrlatt(3,3)
  integer :: kptrlatt_fine(3,3)
  integer :: qrefine(3)

! Real(dp)
  real(dp) :: bmass
  real(dp) :: conf_power_fact_disp
  real(dp) :: conf_power_fact_strain
  real(dp) :: delta_df
  real(dp) :: energy_reference
  real(dp) :: fit_cutoff
  real(dp) :: temperature
  real(dp) :: rifcsph
  real(dp) :: conf
  real(dp) :: acell(3)
  real(dp) :: strain(6)
  real(dp) :: strtarget(6)
  real(dp) :: conf_cutoff_strain(6)
  real(dp) :: rprim(3,3)
  real(dp) :: q1shft(3,4)

! Integer arrays
  integer, allocatable :: atifc(:)
  ! atifc(natom)
  integer, allocatable :: fit_fixcoeff(:)
  ! fit_fixcoeffs(fit_nfixcoeff)

  integer, allocatable :: qmass(:)
  ! qmass(nnos)

! Real arrays
  real(dp), allocatable :: coefficients(:) 
  ! coefficients(ncoeff)

  real(dp), allocatable :: conf_cutoff_disp(:) 
  ! conf_cuttoff(natom)

  real(dp), allocatable :: qnrml1(:) 
  ! qnrml1(nph1l)

  real(dp), allocatable :: qnrml2(:)
  ! qnrml1(nph1l)

  real(dp), allocatable :: qph1l(:,:) 
  ! qph1l(3,nph1l)

  real(dp), allocatable :: qph2l(:,:)
  ! qph2l(3,nph2l)

 end type multibinit_dataset_type
!!***

contains 
!!***

!!****f* m_multibinit_dataset/multibinit_dtset_free
!!
!! NAME
!!  multibinit_dtset_free
!!
!! FUNCTION
!!  deallocate remaining arrays in the multibinit_dtset datastructure
!!
!! INPUTS
!!  multibinit_dtset <type(multibinit_dataset_type)> = multibinit_dataset structure
!!
!! OUTPUTS
!!  multibinit_dtset <type(multibinit_dataset_type)> = multibinit_dataset structure
!!
!! PARENTS
!!      multibinit
!!
!! CHILDREN
!!
!! SOURCE

subroutine multibinit_dtset_free(multibinit_dtset)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'multibinit_dtset_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(multibinit_dataset_type), intent(inout) :: multibinit_dtset

! *************************************************************************
 
 if (allocated(multibinit_dtset%atifc))  then
   ABI_DEALLOCATE(multibinit_dtset%atifc)
 end if
 if (allocated(multibinit_dtset%conf_cutoff_disp))  then
   ABI_DEALLOCATE(multibinit_dtset%conf_cutoff_disp)
 end if
 if (allocated(multibinit_dtset%fit_fixcoeff))  then
   ABI_DEALLOCATE(multibinit_dtset%fit_fixcoeff)
 end if
 if (allocated(multibinit_dtset%qmass))  then
   ABI_DEALLOCATE(multibinit_dtset%qmass)
 end if
 if (allocated(multibinit_dtset%coefficients))  then
   ABI_DEALLOCATE(multibinit_dtset%coefficients)
 end if
 if (allocated(multibinit_dtset%qnrml1))  then
   ABI_DEALLOCATE(multibinit_dtset%qnrml1)
 end if
 if (allocated(multibinit_dtset%qnrml2))  then
   ABI_DEALLOCATE(multibinit_dtset%qnrml2)
 end if
 if (allocated(multibinit_dtset%qph1l))  then
   ABI_DEALLOCATE(multibinit_dtset%qph1l)
 end if
 if (allocated(multibinit_dtset%qph2l))  then
   ABI_DEALLOCATE(multibinit_dtset%qph2l)
 end if

end subroutine multibinit_dtset_free
!!***

!----------------------------------------------------------------------

!!****f* m_multibinit_dataset/invars10
!!
!! NAME
!! invars10
!!
!! FUNCTION
!! Open input file for the multibinit code, then reads or echoes the input information.
!!
!! INPUTS
!! lenstr=actual length of string
!! natom=number of atoms, needed for atifc
!! string*(*)=string of characters containing all input variables and data
!!
!! OUTPUT
!! multibinit_dtset <type(multibinit_dataset_type)> = datatype with all the input variables
!!
!! NOTES
!! Should be executed by one processor only.
!!
!! PARENTS
!!      multibinit
!!
!! CHILDREN
!!
!! SOURCE

subroutine invars10(multibinit_dtset,lenstr,natom,string)


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
 type(multibinit_dataset_type),intent(inout) :: multibinit_dtset 

!Local variables -------------------------
!Dummy arguments for subroutine 'intagm' to parse input file
!Set routine version number here:
!scalars
 integer,parameter :: vrsddb=100401
 integer :: iatifc,ii,iph1,iph2,jdtset,jj,marr,tread
 character(len=500) :: message
!arrays
 integer,allocatable :: intarr(:)
 real(dp),allocatable :: dprarr(:),work(:)

!*********************************************************************
 marr=3
 ABI_ALLOCATE(intarr,(marr))
 ABI_ALLOCATE(dprarr,(marr))

 jdtset=1

!copy natom to multibinit_dtset
 multibinit_dtset%natom=natom

!=====================================================================
!start reading in dimensions and non-dependent variables
!=====================================================================

!A
 multibinit_dtset%asr=2
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'asr',tread,'INT')
 if(tread==1) multibinit_dtset%asr=intarr(1)
 if(multibinit_dtset%asr<-2.or.multibinit_dtset%asr>5)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'asr is',multibinit_dtset%asr,', but the only allowed values',ch10,&
&   'are 0, 1, 2, 3, 4, 5, -1 or -2 .',ch10,&
&   'Action: correct asr in your input file.'
!  Note : negative values are allowed when the acoustic sum rule
!  is to be applied after the analysis of IFCs
!  3,4 are for rotational invariance (under development)
!  5 is for hermitian imposition of the ASR
   MSG_ERROR(message)
 end if

!B
 multibinit_dtset%brav=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'brav',tread,'INT')
 if(tread==1) multibinit_dtset%brav=intarr(1)
 if(multibinit_dtset%brav/=1)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'brav is',multibinit_dtset%brav,', but the only allowed values',ch10,&
&   'are 1 for multibinit (not implemented) .',ch10,&
&   'Action: correct brav in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%bmass=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'bmass',tread,'DPR')
 if(tread==1) multibinit_dtset%bmass=dprarr(1)
 if(multibinit_dtset%bmass<0)then
   write(message, '(a,f10.2,a,a,a,a,a)' )&
&   'bmass is',multibinit_dtset%bmass,', but the only allowed values',ch10,&
&   'is superior to 0.',ch10,&
&   'Action: correct bmass in your input file.'
   MSG_ERROR(message)
 end if


!C
 multibinit_dtset%chneut=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'chneut',tread,'INT')
 if(tread==1) multibinit_dtset%chneut=intarr(1)
 if(multibinit_dtset%chneut<0.or.multibinit_dtset%chneut>2)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'chneut is',multibinit_dtset%chneut,', but the only allowed values',ch10,&
&   'are 0, 1 or 2 .',ch10,&
&   'Action: correct chneut in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%confinement=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'confinement',tread,'INT')
 if(tread==1) multibinit_dtset%confinement=intarr(1)
 if(multibinit_dtset%confinement<0.or.multibinit_dtset%confinement>2)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'confinement is',multibinit_dtset%confinement,', but the only allowed values',ch10,&
&   'are 0, 1 or 2 .',ch10,&
&   'Action: correct confinement in your input file.'
   MSG_ERROR(message)
 end if
 
 multibinit_dtset%conf_power_disp=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'conf_power_disp',tread,'INT')
 if(tread==1) multibinit_dtset%conf_power_disp=intarr(1)
 if(multibinit_dtset%conf_power_disp<0)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'conf_power_disp is',multibinit_dtset%conf_power_disp,', but the only allowed values',ch10,&
&   'positive .',ch10,&
&   'Action: correct conf_power_disp in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%conf_power_strain=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'conf_power_strain',tread,'INT')
 if(tread==1) multibinit_dtset%conf_power_strain=intarr(1)
 if(multibinit_dtset%conf_power_strain<0)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'conf_power_strain is',multibinit_dtset%conf_power_strain,', but the only allowed values',ch10,&
&   'are positive .',ch10,&
&   'Action: correct conf_power_strain in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%conf_power_fact_disp=100
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'conf_power_fact_disp',tread,'DPR')
 if(tread==1) multibinit_dtset%conf_power_fact_disp=dprarr(1)

 multibinit_dtset%conf_power_fact_strain=100
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'conf_power_fact_strain',tread,'DPR')
 if(tread==1) multibinit_dtset%conf_power_fact_strain=dprarr(1)

!D
 multibinit_dtset%dipdip=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dipdip',tread,'INT')
 if(tread==1) multibinit_dtset%dipdip=intarr(1)
 if(multibinit_dtset%dipdip>1.or.multibinit_dtset%dipdip<0)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'dipdip is',multibinit_dtset%dipdip,', but the only allowed values',ch10,&
&   'is 1.',ch10,&
&   'Action: correct dipdip in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%dtion=100
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dtion',tread,'INT')
 if(tread==1) multibinit_dtset%dtion=intarr(1)
 if(multibinit_dtset%dtion<1)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'dtion is',multibinit_dtset%dtion,', but the only allowed values',ch10,&
&   'is superior to 1.',ch10,&
&   'Action: correct dtion in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%delta_df= 1d-02
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'delta_df',tread,'DPR')
 if(tread==1) multibinit_dtset%delta_df=dprarr(1)
 if(multibinit_dtset%delta_df<0)then
   write(message, '(a,es10.2,a,a,a,a,a)' )&
&   'delta_df is',multibinit_dtset%delta_df,', but the only allowed values',ch10,&
&   'are superior to 0  .',ch10,&
&   'Action: correct delta_df in your input file.'
   MSG_ERROR(message)
 end if

!E
 multibinit_dtset%energy_reference= zero
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'energy_reference',tread,'DPR')
 if(tread==1) multibinit_dtset%energy_reference=dprarr(1)

 multibinit_dtset%enunit=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'enunit',tread,'INT')
 if(tread==1) multibinit_dtset%enunit=intarr(1)
 if(multibinit_dtset%enunit<0.or.multibinit_dtset%enunit>2)then
   write(message, '(a,i0,a,a,a,a,a)' )&
&   'enunit is',multibinit_dtset%enunit,', but the only allowed values',ch10,&
&   'are 0, 1 or 2.',ch10,&
&   'Action: correct enunit in your input file.'
   MSG_ERROR(message)
 end if

!F
 multibinit_dtset%fit_option=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'fit_option',tread,'INT')
 if(tread==1) multibinit_dtset%fit_option=intarr(1)
! No test for now
! if(multibinit_dtset%fit_option=0.or.multibinit_dtset%fit_option/=1)then
!   write(message, '(a,i8,a,a,a,a,a)' )&
!&   'fit_option is',multibinit_dtset%fit_option,', but the only allowed values',ch10,&
!&   'are 0 or 1 for multibinit.',ch10,&
!&   'Action: correct fit_option in your input file.'
!   MSG_ERROR(message)
! end if


 multibinit_dtset%fit_ncycle=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'fit_ncycle',tread,'INT')
 if(tread==1) multibinit_dtset%fit_ncycle=intarr(1)
 if(multibinit_dtset%fit_ncycle<0)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'fit_ncycle is',multibinit_dtset%fit_ncycle,', but the only allowed values',ch10,&
&   'are positives for multibinit.',ch10,&
&   'Action: correct fit_ncycle in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%fit_nfixcoeff=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'fit_nfixcoeff',tread,'INT')
 if(tread==1) multibinit_dtset%fit_nfixcoeff=intarr(1)
 if(multibinit_dtset%fit_nfixcoeff<-1)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'fit_nfixcoeff is',multibinit_dtset%fit_nfixcoeff,', but the only allowed values',ch10,&
&   'are -1 or positives for multibinit.',ch10,&
&   'Action: correct fit_nfixcoeff in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%ifcana=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ifcana',tread,'INT')
 if(tread==1) multibinit_dtset%ifcana=intarr(1)
 if(multibinit_dtset%ifcana<0.or.multibinit_dtset%ifcana>1)then
   write(message, '(a,i0,a,a,a,a,a)' )&
&   'ifcana is',multibinit_dtset%ifcana,', but the only allowed values',ch10,&
&   'are 0 or 1.',ch10,&
&   'Action: correct ifcana in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%ifcflag=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ifcflag',tread,'INT')
 if(tread==1) multibinit_dtset%ifcflag=intarr(1)
 if(multibinit_dtset%ifcflag<0.or.multibinit_dtset%ifcflag>1)then
   write(message, '(a,i0,a,a,a,a,a)' )&
&   'ifcflag is',multibinit_dtset%ifcflag,', but the only allowed values',ch10,&
&   'are 0 or 1.',ch10,&
&   'Action: correct ifcflag in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%prtsrlr=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtsrlr',tread,'INT')
 if(tread==1) multibinit_dtset%prtsrlr=intarr(1)
 if(multibinit_dtset%prtsrlr<0.or.multibinit_dtset%prtsrlr>1)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'prtsrlr is',multibinit_dtset%prtsrlr,', but the only allowed values',ch10,&
&   'are 0 or 1.',ch10,&
&   'Action: correct prtsrlr in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%ifcout=2000000 ! or -1 -> max number of ifc
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ifcout',tread,'INT')
 if(tread==1) multibinit_dtset%ifcout=intarr(1)
 if(multibinit_dtset%ifcout<-1)then
   write(message, '(a,i0,a,a,a)' )&
&   'ifcout is',multibinit_dtset%ifcout,', which is lower than -1 (default = all ifc) .',ch10,&
&   'Action: correct ifcout in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%ntime=200
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ntime',tread,'INT')
 if(tread==1) multibinit_dtset%ntime=intarr(1)
 if(multibinit_dtset%ntime<0)then
   write(message, '(a,i0,a,a,a)' )&
&   'ntime is',multibinit_dtset%ntime,', which is lower than 0 .',ch10,&
&   'Action: correct ntime in your input file.'
   MSG_ERROR(message)
 end if


 multibinit_dtset%dynamics=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'dynamics',tread,'INT')
 if(tread==1) multibinit_dtset%dynamics=intarr(1)
 if(multibinit_dtset%dynamics/=0.and.&
&   multibinit_dtset%dynamics/=12.and.multibinit_dtset%dynamics/=13&
&   .and.multibinit_dtset%dynamics/=24.and.multibinit_dtset%dynamics/=25) then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'dynamics is',multibinit_dtset%dynamics,', but the only allowed values',ch10,&
&   'are 12 or  13 (see ionmov in abinit documentation).',ch10,&
&   'Action: correct dynamics in your input file.'
   MSG_ERROR(message)
 end if

!N
 multibinit_dtset%natifc=natom
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'natifc',tread,'INT')
 if(tread==1) multibinit_dtset%natifc=intarr(1)
 if(multibinit_dtset%natifc<0)then
   write(message, '(a,i0,a,a,a)' )&
&   'natifc is',multibinit_dtset%natifc,', which is lower than 0 .',ch10,&
&   'Action: correct natifc in your input file.'
   MSG_ERROR(message)
 end if

 if(multibinit_dtset%natifc>natom)then
   write(message, '(a,i0,a,a,a,i0,a,a,a)' )&
&   'The number of atom ifc in the input files',multibinit_dtset%natifc,',',ch10,&
&   'is larger than the number of atoms',natom,'.',ch10,&
&   'Action: change natifc in the input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%ncoeff=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ncoeff',tread,'INT')
 if(tread==1) multibinit_dtset%ncoeff=intarr(1)
 if(multibinit_dtset%ncoeff<0)then
   write(message, '(a,i0,a,a,a)' )&
&   'ncoeff is',multibinit_dtset%ncoeff,', which is lower than 0 .',ch10,&
&   'Action: correct natifc in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%ng2qpt(:)=0
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'ng2qpt',tread,'INT')
 if(tread==1) multibinit_dtset%ng2qpt(:)=intarr(1:3)
 do ii=1,3
   if(multibinit_dtset%ng2qpt(ii)<0)then
     write(message, '(a,i0,a,i0,a,a,a,i0,a)' )&
&     'ng2qpt(',ii,') is',multibinit_dtset%ng2qpt(ii),', which is lower than 0 .',ch10,&
&     'Action: correct ng2qpt(',ii,') in your input file.'
     MSG_ERROR(message)
   end if
 end do

 multibinit_dtset%n_cell(:)= one
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'n_cell',tread,'INT')
 if(tread==1) multibinit_dtset%n_cell(1:3)=intarr(1:3)
 do ii=1,3
   if(multibinit_dtset%n_cell(ii)<0.or.multibinit_dtset%n_cell(ii)>50)then
     write(message, '(a,i0,a,i0,a,a,a,i0,a)' )&
&     'n_cell(',ii,') is ',multibinit_dtset%n_cell(ii),', which is lower than 0 of superior than 50.',&
&     ch10,'Action: correct n_cell(',ii,') in your input file.'
     MSG_ERROR(message)
   end if
 end do

 multibinit_dtset%ngqpt(:)= one
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'ngqpt',tread,'INT')
 if(tread==1) multibinit_dtset%ngqpt(1:3)=intarr(1:3)
 do ii=1,3
   if(multibinit_dtset%ngqpt(ii)<0)then
     write(message, '(a,i0,a,i0,a,a,a,i0,a)' )&
&     'ngqpt(',ii,') is',multibinit_dtset%ngqpt(ii),', which is lower than 0 .',ch10,&
&     'Action: correct ngqpt(',ii,') in your input file.'
     MSG_ERROR(message)
   end if
 end do

 multibinit_dtset%nph1l=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nph1l',tread,'INT')
 if(tread==1) multibinit_dtset%nph1l=intarr(1)
 if(multibinit_dtset%nph1l<0)then
   write(message, '(a,i0,a,a,a)' )&
&   'nph1l is',multibinit_dtset%nph1l,', which is lower than 0 .',ch10,&
&   'Action: correct nph1l in your input file.'
   MSG_ERROR(message)
 end if
 
 multibinit_dtset%nph2l=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nph2l',tread,'INT')
 if(tread==1) multibinit_dtset%nph2l=intarr(1)
 if(multibinit_dtset%nph2l<0)then
   write(message, '(a,i0,a,a,a)' )&
&   'nph2l is',multibinit_dtset%nph2l,', which is lower than 0 .',ch10,&
&   'Action: correct nph2l in your input file.'
   MSG_ERROR(message)
 end if
 
 multibinit_dtset%nqshft=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nqshft',tread,'INT')
 if(tread==1) multibinit_dtset%nqshft=intarr(1)
 if(multibinit_dtset%nqshft<0 .or. multibinit_dtset%nqshft==3 .or.&
& multibinit_dtset%nqshft>=5 )then
   write(message, '(a,i0,a,a,a,a,a)' )&
&   'nqshft is',multibinit_dtset%nqshft,', but the only allowed values',ch10,&
&   'are 1, 2 or 4 .',ch10,&
&   'Action: correct nqshft in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%nnos=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nnos',tread,'INT')
 if(tread==1) multibinit_dtset%nnos=intarr(1)
 if(multibinit_dtset%nnos<0)then
   write(message, '(a,i0,a,a,a)' )&
&   'nnos is',multibinit_dtset%nnos,', which is lower than 0',ch10,&
&   'Action: correct nnos in your input file.'
   MSG_ERROR(message)
 end if


 multibinit_dtset%nsphere=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nsphere',tread,'INT')
 if(tread==1) multibinit_dtset%nsphere=intarr(1)
 if(multibinit_dtset%nsphere<0)then
   write(message, '(a,i0,a,a,a)' )&
&   'nsphere is',multibinit_dtset%nsphere,', which is lower than 0',ch10,&
&   'Action: correct nsphere in your input file.'
   MSG_ERROR(message)
 end if

!O
 multibinit_dtset%optcell=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'optcell',tread,'INT')
 if(tread==1) multibinit_dtset%optcell=intarr(1)
 if(multibinit_dtset%optcell<0.or.multibinit_dtset%optcell>2)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'optcell is',multibinit_dtset%prtsrlr,', but the only allowed values',ch10,&
&   'are 0, 1 or 2.',ch10,&
&   'Action: correct optcell in your input file.'
   MSG_ERROR(message)
 end if


!P
 multibinit_dtset%prt_model=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prt_model',tread,'INT')
 if(tread==1) multibinit_dtset%prt_model=intarr(1)
 if(multibinit_dtset%prt_model<0.or.multibinit_dtset%prt_model>4)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'prt_model is',multibinit_dtset%prtsrlr,', but the only allowed values',ch10,&
&   'are 0, 1 or 2.',ch10,&
&   'Action: correct prt_model in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%prt_phfrq=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prt_phfrq',tread,'INT')
 if(tread==1) multibinit_dtset%prt_phfrq=intarr(1)
 if(multibinit_dtset%prt_phfrq<0.or.multibinit_dtset%prt_phfrq>2)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'prt_phfrq is',multibinit_dtset%prtsrlr,', but the only allowed values',ch10,&
&   'are 0, 1 or 2.',ch10,&
&   'Action: correct prt_phfrq in your input file.'
   MSG_ERROR(message)
 end if

!Default is no output of the real space IFC to file
 multibinit_dtset%prt_ifc = 0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prt_ifc',tread,'INT')
 if(tread==1) multibinit_dtset%prt_ifc = intarr(1)
 if(multibinit_dtset%prt_ifc < 0 .or. multibinit_dtset%prt_ifc > 1) then
   write(message, '(a,i0,a,a,a,a,a)' )&
&   'prtf_ifc is',multibinit_dtset%prt_ifc,'. The only allowed values',ch10,&
&   'are 0 (no output) or 1 (AI2PS format)',ch10,  &
&   'Action: correct prt_ifc in your input file.'
   MSG_ERROR(message)
 end if

!Default is no output of the 3rd derivative
 multibinit_dtset%strcpling = -1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'strcpling',tread,'INT')
 if(tread==1) multibinit_dtset%strcpling = intarr(1)
 if(multibinit_dtset%strcpling < -1 .or. multibinit_dtset%strcpling > 2) then
   write(message, '(a,i0,a,a,a,a,a,a,a)' )&
&   'prtf_3rd is ',multibinit_dtset%strcpling,'. The only allowed values',ch10,&
&   'are 0 (no computation), 1 (only computation)',ch10,&
&   'or 2 (computation and print in xml file)',ch10,  &
&   'Action: correct strcpling in your input file.'
   MSG_ERROR(message)
 end if

!Q
 multibinit_dtset%qrefine=1 ! default is no refinement
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'qrefine',tread,'INT')
 if(tread==1) multibinit_dtset%qrefine = intarr(1:3)
 do ii=1,3
   if(multibinit_dtset%qrefine(ii) < 1) then
     write(message, '(a,3i0,a,a,a,a,a)' )&
&     'qrefine is',multibinit_dtset%qrefine,' The only allowed values',ch10,&
&     'are integers >= 1 giving the refinement of the ngqpt grid',ch10,&
&     'Action: correct qrefine in your input file.'
     MSG_ERROR(message)
   end if
 end do

!R
 multibinit_dtset%restarxf=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'restarxf',tread,'INT')
 if(tread==1) multibinit_dtset%restarxf=intarr(1)
 if(multibinit_dtset%restarxf < -3 .or. multibinit_dtset%restarxf > 0)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'restarxf is',multibinit_dtset%restarxf,', but the only allowed values',ch10,&
&   'is -2 or 0.',ch10,&
&   'Action: correct restarxf in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%rfmeth=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rfmeth',tread,'INT')
 if(tread==1) multibinit_dtset%rfmeth=intarr(1)
 if(multibinit_dtset%rfmeth<1.or.multibinit_dtset%rfmeth>2)then
   write(message, '(a,i0,a,a,a,a,a)' )&
&   'rfmeth is',multibinit_dtset%rfmeth,', but the only allowed values',ch10,&
&   'are 1 or 2 . ',ch10,&
&   'Action: correct rfmeth in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%rifcsph=zero
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'rifcsph',tread,'DPR')
 if(tread==1) multibinit_dtset%rifcsph=dprarr(1)
 if(multibinit_dtset%rifcsph<-tol12)then
   write(message, '(a,f10.3,a,a,a)' )&
&   'rifcsph is',multibinit_dtset%rifcsph,', which is lower than zero.',ch10,&
&   'Action: correct rifcsph in your input file.'
   MSG_ERROR(message)
 end if

!S
 multibinit_dtset%symdynmat=1
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'symdynmat',tread,'INT')
 if(tread==1) multibinit_dtset%symdynmat=intarr(1)
 if(multibinit_dtset%symdynmat/=0.and.multibinit_dtset%symdynmat/=1)then
   write(message, '(a,i0,a,a,a,a,a)' )&
&   'symdynmat is',multibinit_dtset%symdynmat,'. The only allowed values',ch10,&
&   'are 0, or 1.',ch10,&
&   'Action: correct symdynmat in your input file.'
   MSG_ERROR(message)
 end if

!T
 multibinit_dtset%temperature=325
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'temperature',tread,'DPR')
 if(tread==1) multibinit_dtset%temperature=dprarr(1)
 if(multibinit_dtset%temperature<=0)then
   write(message, '(a,f10.1,a,a,a,a,a)' )&
&   'Temperature is ',multibinit_dtset%temperature,'. The only allowed values',ch10,&
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
 multibinit_dtset%acell= one
 if(3>marr)then
   marr=3
   ABI_DEALLOCATE(intarr)
   ABI_DEALLOCATE(dprarr)
   ABI_ALLOCATE(intarr,(marr))
   ABI_ALLOCATE(dprarr,(marr))
 end if
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'acell',tread,'DPR')
 if(tread==1) multibinit_dtset%acell(1:3)= dprarr(1:3)
 if(any(multibinit_dtset%acell<=tol10))then
    write(message, '(3a)' )&
&       'There is negative on zero value for cell ',ch10,&
&       'Action: change acell in your input file.'
      MSG_ERROR(message) 
 end if

 if(6>marr)then
   marr=6
   ABI_DEALLOCATE(intarr)
   ABI_DEALLOCATE(dprarr)
   ABI_ALLOCATE(intarr,(marr))
   ABI_ALLOCATE(dprarr,(marr))
 end if
 multibinit_dtset%strtarget(1:6) = zero
 call intagm(dprarr,intarr,jdtset,marr,6,string(1:lenstr),'strtarget',tread,'DPR')
 if(tread==1) multibinit_dtset%strtarget(1:6)=dprarr(1:6)


 ABI_ALLOCATE(multibinit_dtset%atifc,(natom))
 multibinit_dtset%atifc(:)=zero
 if(multibinit_dtset%natifc>=1)then
   if(multibinit_dtset%natifc>marr)then
     marr=multibinit_dtset%natifc
     ABI_DEALLOCATE(intarr)
     ABI_DEALLOCATE(dprarr)
     ABI_ALLOCATE(intarr,(marr))
     ABI_ALLOCATE(dprarr,(marr))
   end if
   call intagm(dprarr,intarr,jdtset,marr,multibinit_dtset%natifc,string(1:lenstr),'atifc',tread,'INT')
   if(tread==1) then 
     multibinit_dtset%atifc(1:multibinit_dtset%natifc)= intarr(1:multibinit_dtset%natifc)
   else ! set to the maximum
     do iatifc=1,multibinit_dtset%natifc
       multibinit_dtset%atifc(iatifc) =  iatifc
     end do
   end if
   ABI_MALLOC(work,(natom))
   work(:)=0

   do iatifc=1,multibinit_dtset%natifc
     if(multibinit_dtset%atifc(iatifc)<=0.or.multibinit_dtset%atifc(iatifc)>natom)then
       write(message, '(a,i0,a,a,a,a,a,i0,a,a,a)' )&
&       'For iatifc=',iatifc,', the number of the atom ifc to be ',ch10,&
&       'analysed is not valid : either negative, ',ch10,&
&       'zero, or larger than natom =',natom,'.',ch10,&
&       'Action: change atifc in your input file.'
       MSG_ERROR(message)
     end if
     work(multibinit_dtset%atifc(iatifc))=1
   end do
   multibinit_dtset%atifc(1:natom)=work(:)
   ABI_FREE(work)
 end if

!B

!C
 ABI_ALLOCATE(multibinit_dtset%coefficients,(multibinit_dtset%ncoeff))
 if (multibinit_dtset%ncoeff/=0)then
   if(multibinit_dtset%ncoeff>marr)then
     marr=multibinit_dtset%ncoeff
     ABI_DEALLOCATE(intarr)
     ABI_DEALLOCATE(dprarr)
     ABI_ALLOCATE(intarr,(marr))
     ABI_ALLOCATE(dprarr,(marr))
   end if
   multibinit_dtset%coefficients(:)=zero
   call intagm(dprarr,intarr,jdtset,marr,multibinit_dtset%ncoeff,&
&              string(1:lenstr),'coefficients',tread,'DPR')
   if(tread==1)then
     do ii=1,multibinit_dtset%ncoeff
       multibinit_dtset%coefficients(ii)=dprarr(ii)
     end do
   end if
 end if

 ABI_ALLOCATE(multibinit_dtset%conf_cutoff_disp,(multibinit_dtset%natom))
 if (multibinit_dtset%natom/=0)then
   if(multibinit_dtset%natom>marr)then
     marr=multibinit_dtset%natom
     ABI_DEALLOCATE(intarr)
     ABI_DEALLOCATE(dprarr)
     ABI_ALLOCATE(intarr,(marr))
     ABI_ALLOCATE(dprarr,(marr))
   end if
   multibinit_dtset%conf_cutoff_disp(:)=zero
   call intagm(dprarr,intarr,jdtset,marr,multibinit_dtset%natom,&
&              string(1:lenstr),'conf_cutoff_disp',tread,'DPR')
   if(tread==1)then
     do ii=1,multibinit_dtset%natom
       multibinit_dtset%conf_cutoff_disp(ii)=dprarr(ii)
     end do
   end if
   if(any(multibinit_dtset%conf_cutoff_disp<zero))then
     write(message, '(3a)' )&
&       'There is negative value for conf_cutoff_disp ',ch10,&
&       'Action: change acell in your input file.'
     MSG_ERROR(message) 
   end if
 end if

 if(6>marr)then
   marr=6
   ABI_DEALLOCATE(intarr)
   ABI_DEALLOCATE(dprarr)
   ABI_ALLOCATE(intarr,(marr))
   ABI_ALLOCATE(dprarr,(marr))
 end if
 multibinit_dtset%conf_cutoff_strain(1:6) = zero
 call intagm(dprarr,intarr,jdtset,marr,6,string(1:lenstr),'conf_cutoff_strain',tread,'DPR')
 if(tread==1) multibinit_dtset%conf_cutoff_strain(1:6)=dprarr(1:6)
 if(any(multibinit_dtset%conf_cutoff_disp<zero))then
   write(message, '(3a)' )&
&     'There is negative value for conf_cutoff_strain ',ch10,&
&     'Action: change acell in your input file.'
   MSG_ERROR(message) 
 end if

!D

!E
 multibinit_dtset%eivec=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'eivec',tread,'INT')
 if(tread==1) multibinit_dtset%eivec=intarr(1)
 if(multibinit_dtset%eivec<0.or.multibinit_dtset%eivec>4)then
   write(message, '(a,i0,a,a,a,a,a)' )&
&   'eivec is',multibinit_dtset%eivec,', but the only allowed values',ch10,&
&   'are 0, 1, 2, 3 or 4.',ch10,&
&   'Action: correct eivec in your input file.'
   MSG_ERROR(message)
 end if

!F
 multibinit_dtset%fit_coeff=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'fit_coeff',tread,'INT')
 if(tread==1) multibinit_dtset%fit_coeff=intarr(1)
 if(multibinit_dtset%fit_coeff<0.and.multibinit_dtset%fit_coeff>1)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'fit_coeff is',multibinit_dtset%fit_coeff,', but the only allowed values',ch10,&
&   'are 0 or 1 for multibinit.',ch10,&
&   'Action: correct fit_coeff in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%fit_cutoff=0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'fit_cutoff',tread,'DPR')
 if(tread==1) multibinit_dtset%fit_cutoff=dprarr(1)
 if(multibinit_dtset%fit_cutoff<0)then
   write(message, '(a,i8,a,a,a,a,a)' )&
&   'fit_cutoff is',multibinit_dtset%fit_cutoff,', but the only allowed values',ch10,&
&   'are positives for multibinit.',ch10,&
&   'Action: correct fit_cutoff in your input file.'
   MSG_ERROR(message)
 end if

 multibinit_dtset%fit_grid(:)= one
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'fit_grid',tread,'INT')
 if(tread==1) multibinit_dtset%fit_grid(1:3)=intarr(1:3)
 do ii=1,3
   if(multibinit_dtset%fit_grid(ii)<0.or.multibinit_dtset%fit_grid(ii)>20)then
     write(message, '(a,i0,a,i0,a,a,a,i0,a)' )&
&     'fit_grid(',ii,') is ',multibinit_dtset%fit_grid(ii),', which is lower',&
&     ' than 0 of superior than 20.',&
&     ch10,'Action: correct fit_grid(',ii,') in your input file.'
     MSG_ERROR(message)
   end if
 end do

 multibinit_dtset%fit_rangePower(:)= (/3,4/)
 call intagm(dprarr,intarr,jdtset,marr,2,string(1:lenstr),'fit_rangePower',tread,'INT')
 if(tread==1) multibinit_dtset%fit_rangePower(1:2)=intarr(1:2)
 do ii=1,2
   if(multibinit_dtset%fit_rangePower(ii)<0.or.multibinit_dtset%fit_rangePower(ii)>20)then
     write(message, '(a,i0,a,i0,a,a,a,i0,a)' )&
&     'fit_rangePower(',ii,') is ',multibinit_dtset%fit_rangePower(ii),', which is lower',&
&     ' than 0 of superior than 20.',&
&     ch10,'Action: correct fit_rangePower(',ii,') in your input file.'
     MSG_ERROR(message)
   end if
 end do

!G

!H

!I


!J

!K

!L

!M

!N
 ABI_ALLOCATE(multibinit_dtset%fit_fixcoeff,(multibinit_dtset%fit_nfixcoeff))
 if (multibinit_dtset%fit_nfixcoeff >0)then
   if(multibinit_dtset%fit_nfixcoeff>marr)then
     marr=multibinit_dtset%fit_nfixcoeff
     ABI_DEALLOCATE(intarr)
     ABI_ALLOCATE(intarr,(marr))
   end if
   multibinit_dtset%fit_fixcoeff(:)=zero
   call intagm(dprarr,intarr,jdtset,marr,multibinit_dtset%fit_nfixcoeff,&
&              string(1:lenstr),'fit_fixcoeff',tread,'INT')
   if(tread==1)then
     do ii=1,multibinit_dtset%fit_nfixcoeff
       multibinit_dtset%fit_fixcoeff(ii)=intarr(ii)
     end do
   end if
 end if

!O

!P

!Q
 ABI_ALLOCATE(multibinit_dtset%qmass,(multibinit_dtset%nnos))
 multibinit_dtset%qmass(:)= zero
 if(multibinit_dtset%nnos>=1)then
   if(multibinit_dtset%nnos>marr)then
     marr=multibinit_dtset%nnos
     ABI_DEALLOCATE(intarr)
     ABI_DEALLOCATE(dprarr)
     ABI_ALLOCATE(intarr,(marr))
     ABI_ALLOCATE(dprarr,(marr))
   end if
   call intagm(dprarr,intarr,jdtset,marr,multibinit_dtset%nnos,string(1:lenstr),'qmass',tread,'DPR')
   if(tread==1) multibinit_dtset%qmass(:)=dprarr(1:multibinit_dtset%nnos)
 end if

 if (multibinit_dtset%nqshft/=0)then
   if(3*multibinit_dtset%nqshft>marr)then
     marr=3*multibinit_dtset%nqshft
     ABI_DEALLOCATE(intarr)
     ABI_DEALLOCATE(dprarr)
     ABI_ALLOCATE(intarr,(marr))
     ABI_ALLOCATE(dprarr,(marr))
   end if
   multibinit_dtset%q1shft(:,:)=zero
   call intagm(dprarr,intarr,jdtset,marr,3*multibinit_dtset%nqshft, string(1:lenstr),'q1shft',tread,'DPR')
   if(tread==1) multibinit_dtset%q1shft(1:3,1:multibinit_dtset%nqshft)=&
&   reshape(dprarr(1:3*multibinit_dtset%nqshft),(/3,multibinit_dtset%nqshft/))
 end if

 ABI_ALLOCATE(multibinit_dtset%qph1l,(3,multibinit_dtset%nph1l))
 ABI_ALLOCATE(multibinit_dtset%qnrml1,(multibinit_dtset%nph1l))
 if (multibinit_dtset%nph1l/=0)then
   if(4*multibinit_dtset%nph1l>marr)then
     marr=4*multibinit_dtset%nph1l
     ABI_DEALLOCATE(intarr)
     ABI_DEALLOCATE(dprarr)
     ABI_ALLOCATE(intarr,(marr))
     ABI_ALLOCATE(dprarr,(marr))
   end if
   multibinit_dtset%qph1l(:,:)=zero
   multibinit_dtset%qnrml1(:)=zero
   call intagm(dprarr,intarr,jdtset,marr,4*multibinit_dtset%nph1l,string(1:lenstr),'qph1l',tread,'DPR')
   if(tread==1)then
     do iph1=1,multibinit_dtset%nph1l
       do ii=1,3
         multibinit_dtset%qph1l(ii,iph1)=dprarr(ii+(iph1-1)*4)
       end do
       multibinit_dtset%qnrml1(iph1)=dprarr(4+(iph1-1)*4)
       if(abs(multibinit_dtset%qnrml1(iph1))<DDB_QTOL)then
         write(message, '(a,a,a,a,a)' )&
&         'The first list of wavevectors ','should not have non-analytical data.',ch10,&
&         'Action: correct the first list',' of wavevectors in the input file.'
         MSG_ERROR(message)
       end if
     end do
   end if
 end if

 ABI_ALLOCATE(multibinit_dtset%qph2l,(3,multibinit_dtset%nph2l))
 ABI_ALLOCATE(multibinit_dtset%qnrml2,(multibinit_dtset%nph2l))
 if (multibinit_dtset%nph2l/=0)then
   if(4*multibinit_dtset%nph2l>marr)then
     marr=4*multibinit_dtset%nph2l
     ABI_DEALLOCATE(intarr)
     ABI_DEALLOCATE(dprarr)
     ABI_ALLOCATE(intarr,(marr))
     ABI_ALLOCATE(dprarr,(marr))
   end if
   multibinit_dtset%qph2l(:,:)=zero
   multibinit_dtset%qnrml2(:)=zero
   call intagm(dprarr,intarr,jdtset,marr,4*multibinit_dtset%nph2l,string(1:lenstr),'qph2l',tread,'DPR')
   if(tread==1)then
     do iph2=1,multibinit_dtset%nph2l
       do ii=1,3
         multibinit_dtset%qph2l(ii,iph2)=dprarr(ii+(iph2-1)*4)
       end do
       multibinit_dtset%qnrml2(iph2)=dprarr(4+(iph2-1)*4)
       if(abs(multibinit_dtset%qnrml2(iph2))>DDB_QTOL)then
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
 multibinit_dtset%rprim(:,:)= zero
 call intagm(dprarr,intarr,jdtset,marr,9,string(1:lenstr),'rprim',tread,'DPR')
 if(tread==1) then
   multibinit_dtset%rprim(1:3,1:3)= reshape(dprarr(1:9),(/3,3/))
! check new rprimd
   if(all(multibinit_dtset%rprim(1,:)==zero).or.&
&    all(multibinit_dtset%rprim(2,:)==zero).or.all(multibinit_dtset%rprim(3,:)==zero)) then
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
 multibinit_dtset%strain(:)= zero
 call intagm(dprarr,intarr,jdtset,marr,6,string(1:lenstr),'strain',tread,'DPR')
 if(tread==1) multibinit_dtset%strain(1:6)= dprarr(1:6)

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

 if(multibinit_dtset%prtsrlr/=0 .and. multibinit_dtset%ifcflag/=1) then
   write(message, '(3a)' )&
&   'ifcflag must be 1 for the SR/LR decomposition of the phonon frequencies',ch10,&
&   'Action: correct ifcflag in your input file.'
   MSG_ERROR(message)
 end if

!FIXME: add check that if freeze_displ /= 0 then you need to be doing ifc and phonon interpolation

 if (multibinit_dtset%ifcflag > 0 .and. sum(abs(multibinit_dtset%ngqpt)) == 0) then
   write(message, '(3a)' )&
&   'if you want interatomic force constant output, multibinit needs ngqpt input variable ',ch10,&
&   'Action: set ngqpt in your input file.'
   MSG_ERROR(message)
 end if

!check that q-grid refinement is a divisor of ngqpt in each direction
 if(any(multibinit_dtset%qrefine(:) > 1) .and. &
&    any(abs(dmod(dble(multibinit_dtset%ngqpt(1:3))/dble(multibinit_dtset%qrefine(:)),one)) > tol10)) then
   write(message, '(a,3i0,a,a,a,3i8,a,a)' )&
&   'qrefine is',multibinit_dtset%qrefine,' The only allowed values',ch10,&
&   'are integers which are divisors of the ngqpt grid', multibinit_dtset%ngqpt,ch10,&
&   'Action: correct qrefine in your input file.'
   MSG_ERROR(message)
 end if

! check new rprimd
 if(all(multibinit_dtset%acell(:) > one).and.all(multibinit_dtset%rprim(:,:)==zero))then
   write(message, '(3a)' )&
&         ' acell is defined but there is no rprim',ch10,&
&         'Action: add rprim input'
   MSG_BUG(message)
 end if


!check the fit_fixcoeff
 do ii=1,multibinit_dtset%fit_nfixcoeff
   do jj=ii+1,multibinit_dtset%fit_nfixcoeff
     if (multibinit_dtset%fit_fixcoeff(ii) == multibinit_dtset%fit_fixcoeff(jj))then     
       write(message, '(a,I0,a,I0,2a)' )&
&           ' There is two similar numbers for fit_fixcoeff: ',multibinit_dtset%fit_fixcoeff(ii),&
&           ' and ', multibinit_dtset%fit_fixcoeff(jj),ch10,&
&            'Action: change fit_fixcoeff'
       MSG_BUG(message)
     end if
   end do
 end do

end subroutine invars10
!!***

!----------------------------------------------------------------------

!!****f* m_multibinit_dataset/outvars_multibinit
!!
!! NAME
!! outvars_multibinit
!!
!! FUNCTION
!! Open input file for the multibinit code, then
!! echoes the input information.
!!
!! INPUTS
!! multibinit_dtset <type(multibinit_dataset_type)> datatype with all the input variables 
!! nunit=unit number for input or output
!!
!! OUTPUT
!!  (only writing)
!!
!! NOTES
!! Should be executed by one processor only.
!!
!! PARENTS
!!      multibinit
!!
!! CHILDREN
!!
!! SOURCE

subroutine outvars_multibinit (multibinit_dtset,nunit)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'outvars_multibinit'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: nunit
 type(multibinit_dataset_type),intent(in) :: multibinit_dtset

!Local variables -------------------------
!Set routine version number here:
!scalars
 integer :: ii,iph1,iph2,iqshft

!*********************************************************************

!Write the heading
 write(nunit,'(a,80a,a)') ch10,('=',ii=1,80),ch10
 write(nunit, '(a,a)' )&
& ' -outvars_multibinit: echo values of input variables ----------------------',ch10

!The flags
 if(multibinit_dtset%ifcflag/=0)then
   write(nunit,'(a)')' Flags :'
   if(multibinit_dtset%ifcflag/=0)write(nunit,'(3x,a9,3i10)')'  ifcflag',multibinit_dtset%ifcflag
   if(multibinit_dtset%prt_model/=0)write(nunit,'(3x,a9,3i10)')'prt_model',multibinit_dtset%prt_model
   if(multibinit_dtset%prt_phfrq/=0)write(nunit,'(3x,a9,3i10)')'prt_phfrq',multibinit_dtset%prt_phfrq
   if(multibinit_dtset%strcpling/=0)write(nunit,'(3x,a9,3i10)')'  strcpling',multibinit_dtset%strcpling
   if(multibinit_dtset%strcpling==2)write(nunit,'(3x,a9,3es8.2)')'delta_df',multibinit_dtset%delta_df
 end if

 if(multibinit_dtset%dynamics/=0)then
   write(nunit,'(a)')' Molecular Dynamics :'
   write(nunit,'(3x,a9,3F10.1)')'     temp',multibinit_dtset%temperature
   write(nunit,'(3x,a9,3I10.1)')'    ntime',multibinit_dtset%ntime
   write(nunit,'(3x,a9,3i10)')  '    ncell',multibinit_dtset%n_cell
   write(nunit,'(3x,a9,3i10)')  '    dtion',multibinit_dtset%dtion
   if (multibinit_dtset%restarxf /= zero ) then
     write(nunit,'(3x,a9,3i10)')  ' restarxf',multibinit_dtset%restarxf
   end if
   if(multibinit_dtset%dynamics==13)then
     write(nunit,'(3x,a9,3i10)')'  optcell',multibinit_dtset%optcell
     write(nunit,'(3x,a9,3F12.1)')'    bmass',multibinit_dtset%bmass
     write(nunit,'(3x,a9,3I10)')'     nnos',multibinit_dtset%nnos
     write(nunit,'(3x,a12)',advance='no')'    qmass  '
     write(nunit,'(3x,15i10)') (multibinit_dtset%qmass(ii),ii=1,multibinit_dtset%nnos)
   end if
 end if

 if(multibinit_dtset%confinement==1)then
   write(nunit,'(a)')' Confinement information :'
   write(nunit,'(1x,a22,I5.1)')'       conf_power_disp',multibinit_dtset%conf_power_disp
   write(nunit,'(1x,a22,I5.1)')'     conf_power_strain',multibinit_dtset%conf_power_strain
   write(nunit,'(1x,a22,3es16.8)')'  conf_power_fact_disp',multibinit_dtset%conf_power_fact_disp
   write(nunit,'(1x,a22,3es16.8)')'conf_power_fact_strain',multibinit_dtset%conf_power_fact_strain
   write(nunit,'(1x,a22)')'     conf_cutoff_disp'
   write(nunit,'(19x,3es16.8)') (multibinit_dtset%conf_cutoff_disp(ii),ii=1,multibinit_dtset%natom)
   write(nunit,'(1x,a22)')'    conf_cutoff_strain'
   write(nunit,'(19x,3es16.8)') (multibinit_dtset%conf_cutoff_strain(ii),ii=1,6)
 end if

 if(multibinit_dtset%fit_coeff/=0)then
   write(nunit,'(a)')' Fit the coefficients :'
   write(nunit,'(3x,a14,I10.1)')'     fit_coeff',multibinit_dtset%fit_coeff
   write(nunit,'(3x,a14,F10.1)')'    fit_cutoff',multibinit_dtset%fit_cutoff
   write(nunit,'(3x,a14,I10.1)')'    fit_option',multibinit_dtset%fit_option
   write(nunit,'(3x,a14,I10.1)')'    fit_ncycle',multibinit_dtset%fit_ncycle
   write(nunit,'(3x,a14,3i10)') '      fit_grid',multibinit_dtset%fit_grid
   write(nunit,'(3x,a14,2i10)') 'fit_rangePower',multibinit_dtset%fit_rangePower
   write(nunit,'(3x,a14,I10)')  ' fit_nfixcoeff',multibinit_dtset%fit_nfixcoeff
   write(nunit,'(3x,a14)',advance='no')' fit_fixcoeff'
   write(nunit,'(4x,9i6)') (multibinit_dtset%fit_fixcoeff(ii),ii=1,multibinit_dtset%fit_nfixcoeff)
 end if

!Write the general information
 if( multibinit_dtset%rfmeth/=1 .or. &
& multibinit_dtset%enunit/=0 .or. &
& multibinit_dtset%eivec/=0 .or. &
& multibinit_dtset%asr/=0 .or. &
& multibinit_dtset%chneut/=0)then
   write(nunit,'(a)')' Miscellaneous information :'
   if(multibinit_dtset%rfmeth/=1)write(nunit,'(3x,a9,3i10)')'   rfmeth',multibinit_dtset%rfmeth
   if(multibinit_dtset%enunit/=0)write(nunit,'(3x,a9,3i10)')'   enunit',multibinit_dtset%enunit
   if(multibinit_dtset%eivec/=0) write(nunit,'(3x,a9,3i10)')'    eivec',multibinit_dtset%eivec
   if(multibinit_dtset%asr/=0)   write(nunit,'(3x,a9,3i10)')'      asr',multibinit_dtset%asr
   if(multibinit_dtset%chneut/=0)write(nunit,'(3x,a9,3i10)')'   chneut',multibinit_dtset%chneut
 end if


!For interatomic force constant information
 if(multibinit_dtset%ifcflag/=0)then
   write(nunit,'(a)')' Interatomic Force Constants Inputs :'
   write(nunit,'(3x,a9,3i10)')'   dipdip',multibinit_dtset%dipdip
   if(multibinit_dtset%nsphere/=0)write(nunit,'(3x,a9,3i10)')'  nsphere',multibinit_dtset%nsphere
   if(abs(multibinit_dtset%rifcsph)>tol10)write(nunit,'(3x,a9,E16.6)')'  nsphere',multibinit_dtset%rifcsph
   write(nunit,'(3x,a9,3i10)')'   ifcana',multibinit_dtset%ifcana
   write(nunit,'(3x,a9,3i10)')'   ifcout',multibinit_dtset%ifcout
   if(multibinit_dtset%natifc>=1)then
     write(nunit,'(3x,a9,3i10)')'   natifc',multibinit_dtset%natifc
     write(nunit,'(3x,a12)',advance='no')'    atifc   '
     write(nunit,'(3x,15i4)') (multibinit_dtset%atifc(ii)*ii,ii=1,multibinit_dtset%natifc)

   end if
   write(nunit,'(a)')' Description of grid 1 :'
   write(nunit,'(3x,a9,3i10)')'     brav',multibinit_dtset%brav
   write(nunit,'(3x,a9,3i10)')'    ngqpt',multibinit_dtset%ngqpt(1:3)
   write(nunit,'(3x,a9,3i10)')'   nqshft',multibinit_dtset%nqshft
   if (multibinit_dtset%nqshft/=0)then
     write(nunit,'(3x,a9)')'   q1shft'
     do iqshft=1,multibinit_dtset%nqshft
       write(nunit,'(19x,4es16.8)') (multibinit_dtset%q1shft(ii,iqshft),ii=1,3)
     end do
   end if
   if (any(multibinit_dtset%qrefine(:) > 1)) then
     write(nunit,'(3x,a9,3i10)')'  qrefine', multibinit_dtset%qrefine
   end if
 end if


!List of vector 1  (reduced coordinates)
 if(multibinit_dtset%nph1l/=0)then
   write(nunit,'(a)')' First list of wavevector (reduced coord.) :'
   write(nunit,'(3x,a9,3i10)')'    nph1l',multibinit_dtset%nph1l
   write(nunit,'(3x,a9)')'    qph1l'
   do iph1=1,multibinit_dtset%nph1l
     write(nunit,'(19x,3es16.8,2x,es11.3)') &
&     (multibinit_dtset%qph1l(ii,iph1),ii=1,3),multibinit_dtset%qnrml1(iph1)
   end do
 end if

!List of vector 2  (cartesian coordinates)
 if(multibinit_dtset%nph2l/=0)then
   write(nunit,'(a)')' Second list of wavevector (cart. coord.) :'
   write(nunit,'(3x,a9,3i10)')'    nph2l',multibinit_dtset%nph2l
   write(nunit,'(3x,a9)')'    qph2l'
   do iph2=1,multibinit_dtset%nph2l
     write(nunit,'(19x,3es16.8,2x,es11.3)') &
&     (multibinit_dtset%qph2l(ii,iph2),ii=1,3),multibinit_dtset%qnrml2(iph2)
   end do
 end if

 write(nunit,'(a,80a,a)') ch10,('=',ii=1,80),ch10

end subroutine outvars_multibinit
!!***

end module m_multibinit_dataset
!!***

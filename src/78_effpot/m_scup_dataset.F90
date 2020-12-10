!!****m* ABINIT/m_scup_dataset
!! NAME
!!  m_scup_dataset
!!
!! FUNCTION
!!  module with the type of the input variables for scale_up 
!!  when initialized this is a subtype of multibinit_dtset_type
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2020 ABINIT group (AM)
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

module m_scup_dataset

 use defs_basis
 use m_abicore
 use m_errors

 use m_parser, only : intagm
 use m_symtk,  only : matr3inv
 use m_bz_mesh

 implicit none

 private

 public :: scup_dtset_type
 public :: scup_dtset_init
 public :: scup_dtset_free
 public :: outvars_scup
 public :: invars10scup
 public :: scup_kpath_new 
 public :: scup_kpath_print
!!***

!!****t* m_scup_dataset/scup_dtset_type
!! NAME
!! scup_dtset_type
!!
!! FUNCTION
!! The scup_dtset_type structured datatype
!! gathers all input variables and options for the SCALE UP part 
!! that is linked with multibinit
!!
!! SOURCE

 type scup_dtset_type 

!Integer 
 integer :: scup_nspeck
 integer :: scup_ndivsm
 integer :: scup_printniter
 integer :: scup_startpulay 
 integer :: scup_maxscfstep
!Logicals 
 logical :: scup_elec_model
 logical :: scup_initorbocc
 logical :: scup_ismagnetic 
 logical :: scup_istddft    
 logical :: scup_printbands
 logical :: scup_printeigv
 logical :: scup_printeltic 
 logical :: scup_printgeom 
 logical :: scup_printorbocc
 logical(1) :: scup_freezden
!Real 
 real*8   :: scup_tcharge 
 real*8   :: scup_scfmixing
 real*8   :: scup_scfthresh 
 real*8   :: scup_smearing
!Integer Array
 integer :: scup_ksamp(3)

!Real Array 
 real(dp),allocatable :: scup_speck(:,:)

!Kpath Type 
 type(kpath_t) :: scup_kpath

 end type scup_dtset_type
!!***
CONTAINS

!!****f* m_scup_dataset/scup_dtset_init
!!
!! NAME
!! scup_dtset_init
!!
!! FUNCTION
!! Init the scup_dtset type
!!
!! INPUTS
!! 
!!
!! OUTPUT
!! scup_dtset <type(scup_dtset_type)> = datatype with all the input variables
!!
!! NOTES
!! Should be executed by one processor only.
!!
!! PARENTS
!!      m_scup_dataset
!!
!! CHILDREN
!!      scup_kpath%print
!!
!! SOURCE

subroutine scup_dtset_init(scup_dtset)

implicit none

!Arguments -------------------------------
!scalars
 type(scup_dtset_type),intent(inout) :: scup_dtset
!Local variables -------------------------
!scalars
!arrays
!-----------------------------------------

scup_dtset%scup_ndivsm      =  0
scup_dtset%scup_nspeck      =  0 
scup_dtset%scup_elec_model  = .FALSE.
scup_dtset%scup_ksamp       =  (/ 1, 1, 1 /)
scup_dtset%scup_tcharge     =  0
scup_dtset%scup_initorbocc  = .FALSE. 
scup_dtset%scup_ismagnetic  = .FALSE. 
scup_dtset%scup_istddft     = .FALSE.
scup_dtset%scup_printbands  = .FALSE.  
scup_dtset%scup_printeigv   = .FALSE. 
scup_dtset%scup_printeltic  = .FALSE.  
scup_dtset%scup_printgeom   = .FALSE.  
scup_dtset%scup_printniter  =  0  
scup_dtset%scup_printorbocc = .FALSE. 
scup_dtset%scup_freezden    = .FALSE. 
scup_dtset%scup_scfmixing   =  0.3 
scup_dtset%scup_scfthresh   =  tol6
scup_dtset%scup_smearing    =  0.00091873313 ! Room Temperature in Hartree
scup_dtset%scup_startpulay  =  3
scup_dtset%scup_maxscfstep  =  100

end subroutine  scup_dtset_init
!!***

!!****f* m_scup_dataset/scup_dtset_free
!!
!! NAME
!!  scup_dtset_free
!!
!! FUNCTION
!!  deallocate remaining arrays in the scup_dtset datastructure
!!
!! INPUTS
!!  scup_dtset <type(scup_dtset_type)> = scup_dataset structure
!!
!! OUTPUTS
!!  scup_dtset <type(scup_dtset_type)> = scup_dataset structure
!!
!! PARENTS
!!      m_multibinit_dataset
!!
!! CHILDREN
!!      scup_kpath%print
!!
!! SOURCE

subroutine scup_dtset_free(scup_dtset)

 implicit none

!Arguments ------------------------------------
!scalars
 type(scup_dtset_type), intent(inout) :: scup_dtset

! *************************************************************************


if(allocated(scup_dtset%scup_speck))then 
        ABI_FREE(scup_dtset%scup_speck)
endif

call scup_dtset%scup_kpath%free


end subroutine scup_dtset_free
!!***

!!****f* m_scup_dataset/outvars_scup
!!
!! NAME
!! outvars_scup
!!
!! FUNCTION
!! Takes as an input the input dtset for scup and echoes it to 
!! the output
!!
!! INPUTS
!! multibinit_dtset <type(multibinit_dtset_type)> datatype with all the input variables
!! nunit=unit number for input or output
!!
!! OUTPUT
!!  (only writing)
!!
!! NOTES
!! Should be executed by one processor only.
!!
!! PARENTS
!!      m_multibinit_dataset
!!
!! CHILDREN
!!      scup_kpath%print
!!
!! SOURCE

subroutine outvars_scup(scup_dtset,nunit)

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: nunit
 type(scup_dtset_type),intent(in) :: scup_dtset

!Local variables -------------------------
!Set routine version number here:
!scalars
!integers for printing
 integer :: int_inorb=0,int_mgn=0,int_tddft=0,int_pband=0
 integer :: int_peigv=0,int_peltic=0,int_pgeom=0,int_porbocc=0
 integer :: int_freezden=0
!Character for defining fromat string 
!*********************************************************************
   
   !Check if logicals are true, if yes set integers to one for printing
   if(scup_dtset%scup_initorbocc)   int_inorb    =1
   if(scup_dtset%scup_ismagnetic)   int_mgn      =1
   if(scup_dtset%scup_istddft)      int_tddft    =1
   if(scup_dtset%scup_printbands)   int_pband    =1
   if(scup_dtset%scup_printeigv)    int_peigv    =1
   if(scup_dtset%scup_printeltic)   int_peltic   =1
   if(scup_dtset%scup_printgeom)    int_pgeom    =1
   if(scup_dtset%scup_printorbocc)  int_porbocc  =1
   if(scup_dtset%scup_freezden)     int_freezden =1
    
   !Print  
   write(nunit,'(a)')'Variables for SCALE-UP electronic model :'
   write(nunit,'(1x,a16,3I3)')    '      scup_ksamp',scup_dtset%scup_ksamp
   write(nunit,'(1x,a16,F7.3)')   '    scup_tcharge',scup_dtset%scup_tcharge
   write(nunit,'(1x,a16,I3)')     ' scup_initorbocc',int_inorb  
   write(nunit,'(1x,a16,I3)')     ' scup_ismagnetic',int_mgn    
   write(nunit,'(1x,a16,I3)')     '    scup_istddft',int_tddft  
   write(nunit,'(1x,a16,I3)')     ' scup_printbands',int_pband     
   write(nunit,'(1x,a16,I3)')     '  scup_printeigv',int_peigv   
   write(nunit,'(1x,a16,I3)')     ' scup_printeltic',int_peltic   
   write(nunit,'(1x,a16,I3)')     '  scup_printgeom',int_pgeom    
   write(nunit,'(1x,a16,I3)')     ' scup_printniter',scup_dtset%scup_printniter
   write(nunit,'(1x,a16,I3)')     'scup_printorbocc',int_porbocc 
   write(nunit,'(1x,a16,I3)')     '   scup_freezden',int_freezden
   write(nunit,'(1x,a16,I3)')     '     scup_nspeck',scup_dtset%scup_nspeck
   write(nunit,'(1x,a16,I3)')     '     scup_ndivsm',scup_dtset%scup_ndivsm
   write(nunit,'(1x,a16,F7.3)')   '  scup_scfmixing',scup_dtset%scup_scfmixing 
   write(nunit,'(1x,a16,ES10.2)') '  scup_scfthresh',scup_dtset%scup_scfthresh 
   write(nunit,'(1x,a16,ES10.2)') '   scup_smearing',scup_dtset%scup_smearing 
   write(nunit,'(1x,a16,I3)')     ' scup_startpulay',scup_dtset%scup_startpulay
   write(nunit,'(1x,a16,I3)')     ' scup_maxscfstep',scup_dtset%scup_maxscfstep


end subroutine outvars_scup
!!***

!!****f* m_scup_dataset/invars10scup
!!
!! NAME
!! invars10scup
!!
!! FUNCTION
!! Open input file for the multibinit code, then reads or echoes the input information 
!! for SCALE UP part that is linked with multibinit.
!!
!! INPUTS
!! lenstr=actual length of string
!! natom=number of atoms, needed for atifc
!! string*(*)=string of characters containing all input variables and data
!!
!! OUTPUT
!! scup_dtset <type(scup_dtset_type)> = datatype with all the input variables
!!
!! NOTES
!! Should be executed by one processor only.
!!
!! PARENTS
!!      m_multibinit_dataset
!!
!! CHILDREN
!!      scup_kpath%print
!!
!! SOURCE

subroutine invars10scup(scup_dtset,lenstr,string)

implicit none 
!Arguments -------------------------------
!scalars
 integer,intent(in) :: lenstr
 character(len=*),intent(in) :: string
 type(scup_dtset_type),intent(inout) :: scup_dtset

!Local variables -------------------------
!Dummy arguments for subroutine 'intagm' to parse input file
!Set routine version number here:
!scalars
 integer :: ii,jdtset,marr,tread
 character(len=500) :: message
!arrays
 integer,allocatable :: intarr(:)
 real(dp),allocatable :: dprarr(:)
 !tmp integer to transfer to logicals 
 integer :: tmp_int

!*********************************************************************

 marr=30
 ABI_MALLOC(intarr,(marr))
 ABI_MALLOC(dprarr,(marr))

 jdtset=1

!=====================================================================
! Initialize Dataset with default values
!=====================================================================

call scup_dtset_init(scup_dtset) 
tmp_int=0 

!=====================================================================
!start reading in dimensions and non-dependent variables
!=====================================================================

!A 
!B 
!C 
!D 
!E 
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'scup_elec_model',tread,'INT')
 if(tread==1) tmp_int=intarr(1)
 if(tmp_int<0 .or. tmp_int>1 )then
   write(message, '(a,I3,a,a,a,a,a)' )&
&   'scup_elec_model is',tmp_int,', but the only allowed values',ch10,&
&   'are 0 and 1.',ch10,&
&   'Action: correct scup_elec_model in your input file.'
   ABI_ERROR(message)
 end if
 if(tmp_int == 1) scup_dtset%scup_elec_model = .TRUE.
 tmp_int = 0 

!F 
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'scup_freezden',tread,'INT')
 if(tread==1) tmp_int=intarr(1)
 if(tmp_int<0 .or. tmp_int>1 )then
   write(message, '(a,I3,a,a,a,a,a)' )&
&   'scup_freezden is',tmp_int,', but the only allowed values',ch10,&
&   'are 0 and 1.',ch10,&
&   'Action: correct scup_freezden in your input file.'
   ABI_ERROR(message)
 end if
 if(tmp_int == 1) scup_dtset%scup_freezden = .TRUE.
 tmp_int = 0 

 !G 
!H 
!I 
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'scup_initorbocc',tread,'INT')
 if(tread==1) tmp_int=intarr(1)
 if(tmp_int<0 .or. tmp_int>1 )then
   write(message, '(a,I3,a,a,a,a,a)' )&
&   'scup_initorbocc is',tmp_int,', but the only allowed values',ch10,&
&   'are 0 and 1.',ch10,&
&   'Action: correct scup_initorbocc in your input file.'
   ABI_ERROR(message)
 end if
 if(tmp_int == 1) scup_dtset%scup_initorbocc = .TRUE.
 tmp_int = 0 

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'scup_ismagnetic',tread,'INT')
 if(tread==1) tmp_int=intarr(1)
 if(tmp_int<0 .or. tmp_int>1 )then
   write(message, '(a,I3,a,a,a,a,a)' )&
&   'scup_ismagnetic is',tmp_int,', but the only allowed values',ch10,&
&   'are 0 and 1.',ch10,&
&   'Action: correct scup_ismagnetic in your input file.'
   ABI_ERROR(message)
 end if
 if(tmp_int == 1) scup_dtset%scup_ismagnetic = .TRUE.
 tmp_int = 0 

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'scup_istddft',tread,'INT')
 if(tread==1) tmp_int=intarr(1)
 if(tmp_int<0 .or. tmp_int>1 )then
   write(message, '(a,I3,a,a,a,a,a)' )&
&   'scup_istddft is',tmp_int,', but the only allowed values',ch10,&
&   'are 0 and 1.',ch10,&
&   'Action: correct scup_istddft in your input file.'
   ABI_ERROR(message)
 end if
 if(tmp_int == 1) scup_dtset%scup_istddft = .TRUE.
 tmp_int = 0 


!J
!K 
 call intagm(dprarr,intarr,jdtset,marr,3,string(1:lenstr),'scup_ksamp',tread,'INT')
 if(tread==1) scup_dtset%scup_ksamp(1:3)=intarr(1:3)
 do ii=1,3
   if(scup_dtset%scup_ksamp(ii)<1)then
     write(message, '(a,i0,a,i0,a,a,a,i0,a)' )&
&     'scup_ksamp(',ii,') is',scup_dtset%scup_ksamp(ii),', which is lower than 1 .',ch10,&
&     'Action: correct scup_ksamp(',ii,') in your input file.'
     ABI_ERROR(message)
   end if
 end do

!L 
!M
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'scup_maxscfstep',tread,'INT')
 if(tread==1) scup_dtset%scup_maxscfstep=intarr(1)
 if(scup_dtset%scup_maxscfstep<=0)then
   write(message, '(a,I3,a,a,a,a,a)' )&
&   'scup_maxscfstep is',scup_dtset%scup_maxscfstep,', but the only allowed values',ch10,&
&   'greater than 0',ch10,&
&   'Action: correct scup_maxscfstep in your input file.'
   ABI_ERROR(message)
 end if

!N

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'scup_ndivsm',tread,'INT')
 if(tread==1) scup_dtset%scup_ndivsm=intarr(1)
 if(scup_dtset%scup_ndivsm<0 )then
   write(message, '(a,I3,a,a,a,a,a)' )&
&   'scup_ndivsm is',scup_dtset%scup_ndivsm,', but the only allowed values',ch10,&
&   'are positiv.',ch10,&
&   'Action: correct scup_ndivsm in your input file.'
   ABI_ERROR(message)
 end if


 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'scup_nspeck',tread,'INT')
 if(tread==1) scup_dtset%scup_nspeck=intarr(1)
 if(scup_dtset%scup_nspeck<0 )then
   write(message, '(a,I3,a,a,a,a,a)' )&
&   'scup_nspeck is',scup_dtset%scup_nspeck,', but the only allowed values',ch10,&
&   'are positiv.',ch10,&
&   'Action: correct scup_nspeck in your input file.'
   ABI_ERROR(message)
 end if

!O
!P 
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'scup_printbands',tread,'INT')
 if(tread==1) tmp_int=intarr(1)
 if(tmp_int<0 .or. tmp_int>1 )then
   write(message, '(a,I3,a,a,a,a,a)' )&
&   'scup_printbands is',tmp_int,', but the only allowed values',ch10,&
&   'are 0 and 1.',ch10,&
&   'Action: correct scup_printbands in your input file.'
   ABI_ERROR(message)
 end if
 if(tmp_int == 1) scup_dtset%scup_printbands = .TRUE.
 tmp_int = 0 

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'scup_printeigv',tread,'INT')
 if(tread==1) tmp_int=intarr(1)
 if(tmp_int<0 .or. tmp_int>1 )then
   write(message, '(a,I3,a,a,a,a,a)' )&
&   'scup_printeigv is',tmp_int,', but the only allowed values',ch10,&
&   'are 0 and 1.',ch10,&
&   'Action: correct scup_printeigv in your input file.'
   ABI_ERROR(message)
 end if
 if(tmp_int == 1) scup_dtset%scup_printeigv = .TRUE.
 tmp_int = 0 

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'scup_printeltic',tread,'INT')
 if(tread==1) tmp_int=intarr(1)
 if(tmp_int<0 .or. tmp_int>1 )then
   write(message, '(a,I3,a,a,a,a,a)' )&
&   'scup_printeltic is',tmp_int,', but the only allowed values',ch10,&
&   'are 0 and 1.',ch10,&
&   'Action: correct scup_printeltic in your input file.'
   ABI_ERROR(message)
 end if
 if(tmp_int == 1) scup_dtset%scup_printeltic = .TRUE.
 tmp_int = 0 

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'scup_printgeom',tread,'INT')
 if(tread==1) tmp_int=intarr(1)
 if(tmp_int<0 .or. tmp_int>1 )then
   write(message, '(a,I3,a,a,a,a,a)' )&
&   'scup_printgeom is',tmp_int,', but the only allowed values',ch10,&
&   'are 0 and 1.',ch10,&
&   'Action: correct scup_printgeom in your input file.'
   ABI_ERROR(message)
 end if
 if(tmp_int == 1) scup_dtset%scup_printgeom = .TRUE.
 tmp_int = 0 

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'scup_printniter',tread,'INT')
 if(tread==1) scup_dtset%scup_printniter=intarr(1)
 if(scup_dtset%scup_printniter<0 )then
   write(message, '(a,I3,a,a,a,a,a)' )&
&   'scup_printniter is',scup_dtset%scup_printniter,', but the only allowed values',ch10,&
&   'are positiv',ch10,&
&   'Action: correct scup_printniter in your input file.'
   ABI_ERROR(message)
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'scup_printorbocc',tread,'INT')
 if(tread==1) tmp_int=intarr(1)
 if(tmp_int<0 .or. tmp_int>1 )then
   write(message, '(a,I3,a,a,a,a,a)' )&
&   'scup_printorbocc is',tmp_int,', but the only allowed values',ch10,&
&   'are 0 and 1.',ch10,&
&   'Action: correct scup_printorbocc in your input file.'
   ABI_ERROR(message)
 end if
 if(tmp_int == 1) scup_dtset%scup_printorbocc = .TRUE.
 tmp_int = 0 

!Q
!R
!S

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'scup_startpulay',tread,'INT')
 if(tread==1) scup_dtset%scup_startpulay=intarr(1)
 if(scup_dtset%scup_startpulay<3)then
   write(message, '(a,I3,a,a,a,a,a)' )&
&   'scup_startpulay is',scup_dtset%scup_startpulay,', but the only allowed values',ch10,&
&   'are greater than 3',ch10,&
&   'Action: correct scup_startpulay in your input file.'
   ABI_ERROR(message)
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'scup_scfmixing',tread,'DPR')
 if(tread==1) scup_dtset%scup_scfmixing=dprarr(1)
 if(scup_dtset%scup_scfmixing<0)then
   write(message, '(a,f10.2,a,a,a,a,a)' )&
&   'scup_scfmixing is',scup_dtset%scup_scfmixing,', but the only allowed value',ch10,&
&   'is superior to 0.',ch10,&
&   'Action: correct scup_scfmixing in your input file.'
   ABI_ERROR(message)
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'scup_scfthresh',tread,'DPR')
 if(tread==1) scup_dtset%scup_scfthresh=dprarr(1)
 if(scup_dtset%scup_scfthresh <= 0)then
   write(message, '(a,f10.2,a,a,a,a,a)' )&
&   'scup_scfthresh is',scup_dtset%scup_scfthresh,', but the only allowed value',ch10,&
&   'is superior to 0.',ch10,&
&   'Action: correct scup_scfthresh in your input file.'
   ABI_ERROR(message)
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'scup_smearing',tread,'DPR')
 if(tread==1) scup_dtset%scup_smearing=dprarr(1)
 if(scup_dtset%scup_smearing < 0)then
   write(message, '(a,f10.2,a,a,a,a,a)' )&
&   'scup_smearing is',scup_dtset%scup_smearing,', but the only allowed value',ch10,&
&   'is superior to or equal to 0.',ch10,&
&   'Action: correct scup_smearing in your input file.'
   ABI_ERROR(message)
 end if

!T
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'scup_tcharge',tread,'DPR')
 if(tread==1) scup_dtset%scup_tcharge=dprarr(1)
 if(scup_dtset%scup_tcharge<0)then
   write(message, '(a,f10.2,a,a,a,a,a)' )&
&   'scup_tcharge is',scup_dtset%scup_tcharge,', but the only allowed value',ch10,&
&   'is superior to 0.',ch10,&
&   'Action: correct scup_tcharge in your input file.'
   ABI_ERROR(message)
 end if


!U 
!V 
!W
!X
!Y
!Z 

!=====================================================================
!start reading dimension dependent variables
!=====================================================================


!A
!B
!C
!D
!E
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
!R
!S
   !Allocate
   if(scup_dtset%scup_printbands)then 
     ABI_MALLOC(scup_dtset%scup_speck,(3,scup_dtset%scup_nspeck))
    
   call intagm(dprarr,intarr,jdtset,marr,3*scup_dtset%scup_nspeck,string(1:lenstr),'scup_speck',tread,'DPR')
     if(tread==1)then
       scup_dtset%scup_speck(:,:)=reshape( dprarr(1:3*scup_dtset%scup_nspeck), [3,scup_dtset%scup_nspeck])
     else
       write(message,'(5a)') &
&       'When scup_printbands is asked, scup_speck must be initialized ',ch10,&
&       'in the input file, which is not the case.',ch10,&
&       'Action: initialize scup_speck in your input file, or change printbands.'
       ABI_ERROR(message)
     end if
   end if


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

 ABI_FREE(dprarr)
 ABI_FREE(intarr)

end subroutine invars10scup
!!***

!!****f* m_scup_dataset/scup_kpath_new
!!
!! NAME
!! scup_kpath_init
!!
!! FUNCTION
!! Initialize the kpath and all other variables SCALE UP needs to plot electronic bands 
!! along kpath
!!
!! INPUTS
!!
!! speck = array with special k-points along the path 
!! gprimd = reciprocal latice vectors of cell 
!! ndivsm = number of divisions for smallest segment 
!!
!!
!! OUTPUT
!! scup_kpath <type(kpath_t)> = kpath_t with all information about kpath
!!
!! NOTES
!! Should be executed by one processor only.
!!
!! PARENTS
!!      multibinit
!!
!! CHILDREN
!!  
!!  m_bz_mesh/kpath_new
!!
!! SOURCE

subroutine scup_kpath_new(speck,rprimd,ndivsm,scup_kpath)

 implicit none
!Arguments -------------------------------
!scalars
 integer,intent(in) :: ndivsm

!arrays 
 real(dp),intent(in) :: speck(:,:),rprimd(3,3)
 type(kpath_t), intent(out) :: scup_kpath

!Local variables -------------------------
!Dummy arguments for subroutine 'intagm' to parse input file
!Set routine version number here:
!scalars
 integer :: nspeck 
!arrays
 integer,allocatable :: ndivs_tmp(:)
 real(dp) :: gprimd(3,3)

!*********************************************************************

!Get gprimd 
call matr3inv(rprimd,gprimd)

!Create Kpath 
scup_kpath = kpath_new(speck,gprimd,ndivsm)

!Change size of scup_kpath%ndivs(:) variable 
!from nspeck-1 to nspeck and put 1 to first entry
nspeck = size(speck,2)
ABI_MALLOC(ndivs_tmp,(nspeck))

!First entry is always 1
 ndivs_tmp(1) = 1 
!Copy point/per segments
 ndivs_tmp(2:) = scup_kpath%ndivs(:) 

!Delete original segments 
ABI_FREE(scup_kpath%ndivs)
!Put new path 
ABI_MALLOC(scup_kpath%ndivs,(nspeck))
scup_kpath%ndivs = ndivs_tmp
!Free temporary array
ABI_FREE(ndivs_tmp)

end subroutine scup_kpath_new
!!***

!!****f* m_scup_dataset/scup_kpath_print
!!
!! NAME
!! scup_kpath_print
!!
!! FUNCTION
!! Print info of kpath provide to SCALE UP 
!!
!! INPUTS
!!
!! scup_kpath<tpye(kpath_t) = kpath_t with all information about kpath  
!!
!! OUTPUT
!! Only Printing
!!
!! NOTES
!! Should be executed by one processor only.
!!
!! PARENTS
!!      m_multibinit_driver
!!
!! CHILDREN
!!      scup_kpath%print
!!
!! SOURCE

subroutine scup_kpath_print(scup_kpath)

 implicit none 
!Arguments ------------------------------------
!scalars
!arrays 
 type(kpath_t), intent(in) :: scup_kpath 
!Local variables-------------------------------
 integer :: unt

! *************************************************************************

 unt = std_out 
 
 write(unt,'(a)') ch10
 write(unt,'(4a)') ' scup_printbands = 1. Printing of electronic bands active',ch10,&
&                  ' Kpath information below:',ch10

call scup_kpath%print

end subroutine scup_kpath_print
!!***

end module m_scup_dataset
!!***

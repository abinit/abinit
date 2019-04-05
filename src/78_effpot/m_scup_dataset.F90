!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_scup_dataset
!! NAME
!!  m_scup_dataset
!!
!! FUNCTION
!!  module with the type of the input variables for scale_up 
!!  when initialized this is a subtype of multibinit_dtset_type
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2018 ABINIT group (AM)
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
 use m_ddb,    only : DDB_QTOL
 use m_bz_mesh

 implicit none

 private

 public :: scup_dtset_type
 public :: scup_dtset_init
 public :: scup_dtset_free
 public :: outvars_scup
 public :: invars10scup


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
!Logicals 
 logical*1 :: scup_elec_model
 logical*1 :: scup_initorbocc
 logical*1 :: scup_ismagnetic 
 logical*1 :: scup_istddft    
 logical*1 :: scup_printbands
 logical*1 :: scup_printeigv
 logical*1 :: scup_printeltic 
 logical*1 :: scup_printgeom 
 logical*1 :: scup_printorbocc
!Real 
 real*8   :: scup_tcharge 

!Integer Array
 integer :: scup_ksamp(3)

!Real Array 
 real,allocatable :: scup_speck(:,:)

 end type scup_dtset_type
!!***
contains

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
!!
!! CHILDREN
!!
!! SOURCE

subroutine scup_dtset_init(scup_dtset)

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
scup_dtset%scup_printorbocc = .FALSE. 



end subroutine  scup_dtset_init


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
!!      multibinit
!!
!! CHILDREN
!!
!! SOURCE

subroutine scup_dtset_free(scup_dtset)

 implicit none

!Arguments ------------------------------------
!scalars
 type(scup_dtset_type), intent(inout) :: scup_dtset

! *************************************************************************


if(allocated(scup_dtset%scup_speck))then 
        ABI_DEALLOCATE(scup_dtset%scup_speck)
endif


end subroutine scup_dtset_free


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
!!      multibinit
!!
!! CHILDREN
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
!Character for defining fromat string 
 character(len=32) :: string
!*********************************************************************
   
   !Check if logicals are true, if yes set integers to one for printing
   if(scup_dtset%scup_initorbocc)   int_inorb   =1
   if(scup_dtset%scup_ismagnetic)   int_mgn     =1
   if(scup_dtset%scup_istddft)      int_tddft   =1
   if(scup_dtset%scup_printbands)   int_pband   =1
   if(scup_dtset%scup_printeigv)    int_peigv   =1
   if(scup_dtset%scup_printeltic)   int_peltic  =1
   if(scup_dtset%scup_printgeom)    int_pgeom   =1
   if(scup_dtset%scup_printorbocc)  int_porbocc =1
   
   !Debug format string for writing kpoints to output
   !Write format string for special kpoints 
   !if(allocated(scup_dtset%scup_speck))then 
   !   write (string, '( "(1x,a16 ", I4, "(F6.2))" )' )  scup_dtset%scup_nspeck
   !endif
   
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
   write(nunit,'(1x,a16,I3)')     'scup_printorbocc',int_porbocc 
   write(nunit,'(1x,a16,I3)')     '     scup_nspeck',scup_dtset%scup_nspeck
   write(nunit,'(1x,a16,I3)')     '     scup_ndivsm',scup_dtset%scup_ndivsm
   !Debug write kpoints to output 
   !if(allocated(scup_dtset%scup_speck))then 
   !   write(nunit,string)     '      scup_speck',scup_dtset%scup_speck
   !endif


end subroutine outvars_scup


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
!!      multibinit
!!
!! CHILDREN
!!
!! SOURCE

subroutine invars10scup(scup_dtset,lenstr,string)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: lenstr
 character(len=*),intent(in) :: string
 type(scup_dtset_type),intent(inout) :: scup_dtset

!Local variables -------------------------
!Dummy arguments for subroutine 'intagm' to parse input file
!Set routine version number here:
!scalars
 integer :: ii,jdtset,jj,marr,tread
 character(len=500) :: message
!arrays
 integer,allocatable :: intarr(:)
 real(dp),allocatable :: dprarr(:),work(:)
 !tmp integer to transfer to logicals 
 integer :: tmp_int

!*********************************************************************

 marr=30
 ABI_ALLOCATE(intarr,(marr))
 ABI_ALLOCATE(dprarr,(marr))

 jdtset=1

!=====================================================================
! Initialize Dataset with default values
!=====================================================================

call scup_dtset_init(scup_dtset) 

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
   MSG_ERROR(message)
 end if
 if(tmp_int == 1) scup_dtset%scup_elec_model = .TRUE.
 tmp_int = 0 

!F 
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
   MSG_ERROR(message)
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
   MSG_ERROR(message)
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
   MSG_ERROR(message)
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
     MSG_ERROR(message)
   end if
 end do

!L 
!M 
!N

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'scup_ndivsm',tread,'INT')
 if(tread==1) scup_dtset%scup_ndivsm=intarr(1)
 if(scup_dtset%scup_ndivsm<0 )then
   write(message, '(a,I3,a,a,a,a,a)' )&
&   'scup_ndivsm is',scup_dtset%scup_ndivsm,', but the only allowed values',ch10,&
&   'are positiv.',ch10,&
&   'Action: correct scup_ndivsm in your input file.'
   MSG_ERROR(message)
 end if


 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'scup_nspeck',tread,'INT')
 if(tread==1) scup_dtset%scup_nspeck=intarr(1)
 if(scup_dtset%scup_nspeck<0 )then
   write(message, '(a,I3,a,a,a,a,a)' )&
&   'scup_nspeck is',scup_dtset%scup_nspeck,', but the only allowed values',ch10,&
&   'are positiv.',ch10,&
&   'Action: correct scup_nspeck in your input file.'
   MSG_ERROR(message)
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
   MSG_ERROR(message)
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
   MSG_ERROR(message)
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
   MSG_ERROR(message)
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
   MSG_ERROR(message)
 end if
 if(tmp_int == 1) scup_dtset%scup_printgeom = .TRUE.
 tmp_int = 0 

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'scup_printorbocc',tread,'INT')
 if(tread==1) tmp_int=intarr(1)
 if(tmp_int<0 .or. tmp_int>1 )then
   write(message, '(a,I3,a,a,a,a,a)' )&
&   'scup_printorbocc is',tmp_int,', but the only allowed values',ch10,&
&   'are 0 and 1.',ch10,&
&   'Action: correct scup_printorbocc in your input file.'
   MSG_ERROR(message)
 end if
 if(tmp_int == 1) scup_dtset%scup_printorbocc = .TRUE.
 tmp_int = 0 

!Q
!R
!S 
!T
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'scup_tcharge',tread,'DPR')
 if(tread==1) scup_dtset%scup_tcharge=dprarr(1)
 if(scup_dtset%scup_tcharge<0)then
   write(message, '(a,f10.2,a,a,a,a,a)' )&
&   'scup_tcharge is',scup_dtset%scup_tcharge,', but the only allowed value',ch10,&
&   'is superior to 0.',ch10,&
&   'Action: correct scup_tcharge in your input file.'
   MSG_ERROR(message)
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
     ABI_ALLOCATE(scup_dtset%scup_speck,(3,scup_dtset%scup_nspeck))
    
   call intagm(dprarr,intarr,jdtset,marr,3*scup_dtset%scup_nspeck,string(1:lenstr),'scup_speck',tread,'DPR')
     if(tread==1)then
       scup_dtset%scup_speck(:,:)=reshape( dprarr(1:3*scup_dtset%scup_nspeck), [3,scup_dtset%scup_nspeck])
     else
       write(message,'(5a)') &
&       'When scup_printbands is asked, scup_speck must be initialized ',ch10,&
&       'in the input file, which is not the case.',ch10,&
&       'Action: initialize scup_speck in your input file, or change printbands.'
       MSG_ERROR(message)
     end if
   end if


!T
!U 
!V 
!W
!X
!Y
!Z 


end subroutine invars10scup

end module m_scup_dataset

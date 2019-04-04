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
 integer :: scup_elec_model
 integer :: scup_initorbocc
 integer :: scup_ismagnetic 
 integer :: scup_istddft    
 integer :: scup_printbands
 integer :: scup_printeigv
 integer :: scup_printeltic 
 integer :: scup_printgeom 
 integer :: scup_printorbocc
!Real 
 real*8   :: scup_tcharge 

!Integer Array
 integer :: scup_ksamp(3)

 end type scup_dtset_type
!!***


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

scup_dtset%scup_elec_model  = 0
scup_dtset%scup_ksamp       = (/ 1, 1, 1 /)
scup_dtset%scup_tcharge     = 0 
scup_dtset%scup_initorbocc  = 0 
scup_dtset%scup_ismagnetic  = 0 
scup_dtset%scup_istddft     = 0
scup_dtset%scup_printbands  = 0  
scup_dtset%scup_printeigv   = 0 
scup_dtset%scup_printeltic  = 0  
scup_dtset%scup_printgeom   = 0  
scup_dtset%scup_printorbocc = 0 



end subroutine  scup_dtset_init


!!****f* m_scup_dataset/scup_dtset_free
!!
!! NAME
!!  scup_dtset_free
!!
!! FUNCTION
!!  deallocate remaining arrays in the scup_dtset datastructure
!!  At the moment no allocations in scup_dtset so it's empty
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
 integer :: ii,iph1,iph2,iqshft

!*********************************************************************

   write(nunit,'(a)')'Variables for SCALE-UP electronic model :'
   write(nunit,'(1x,a16,3I3)')    '      scup_ksamp',scup_dtset%scup_ksamp
   write(nunit,'(1x,a16,F7.3)')   '    scup_tcharge',scup_dtset%scup_tcharge
   write(nunit,'(1x,a16,I3)')     ' scup_initorbocc',scup_dtset%scup_initorbocc
   write(nunit,'(1x,a16,I3)')     ' scup_ismagnetic',scup_dtset%scup_ismagnetic
   write(nunit,'(1x,a16,I3)')     '    scup_istddft',scup_dtset%scup_istddft
   write(nunit,'(1x,a16,I3)')     ' scup_printbands',scup_dtset%scup_printbands    
   write(nunit,'(1x,a16,I3)')     '  scup_printeigv',scup_dtset%scup_printeigv   
   write(nunit,'(1x,a16,I3)')     ' scup_printeltic',scup_dtset%scup_printeltic   
   write(nunit,'(1x,a16,I3)')     '  scup_printgeom',scup_dtset%scup_printgeom    
   write(nunit,'(1x,a16,I3)')     'scup_printorbocc',scup_dtset%scup_printorbocc 

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
 if(tread==1) scup_dtset%scup_elec_model=intarr(1)
 if(scup_dtset%scup_elec_model<0 .or. scup_dtset%scup_elec_model>1 )then
   write(message, '(a,I3,a,a,a,a,a)' )&
&   'scup_elec_model is',scup_dtset%scup_elec_model,', but the only allowed values',ch10,&
&   'are 0 and 1.',ch10,&
&   'Action: correct scup_elec_model in your input file.'
   MSG_ERROR(message)
 end if

!F 
!G 
!H 
!I 
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'scup_initorbocc',tread,'INT')
 if(tread==1) scup_dtset%scup_initorbocc=intarr(1)
 if(scup_dtset%scup_initorbocc<0 .or. scup_dtset%scup_initorbocc>1 )then
   write(message, '(a,I3,a,a,a,a,a)' )&
&   'scup_initorbocc is',scup_dtset%scup_initorbocc,', but the only allowed values',ch10,&
&   'are 0 and 1.',ch10,&
&   'Action: correct scup_initorbocc in your input file.'
   MSG_ERROR(message)
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'scup_ismagnetic',tread,'INT')
 if(tread==1) scup_dtset%scup_ismagnetic=intarr(1)
 if(scup_dtset%scup_ismagnetic<0 .or. scup_dtset%scup_ismagnetic>1 )then
   write(message, '(a,I3,a,a,a,a,a)' )&
&   'scup_ismagnetic is',scup_dtset%scup_ismagnetic,', but the only allowed values',ch10,&
&   'are 0 and 1.',ch10,&
&   'Action: correct scup_ismagnetic in your input file.'
   MSG_ERROR(message)
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'scup_istddft',tread,'INT')
 if(tread==1) scup_dtset%scup_istddft=intarr(1)
 if(scup_dtset%scup_istddft<0 .or. scup_dtset%scup_istddft>1 )then
   write(message, '(a,I3,a,a,a,a,a)' )&
&   'scup_istddft is',scup_dtset%scup_istddft,', but the only allowed values',ch10,&
&   'are 0 and 1.',ch10,&
&   'Action: correct scup_istddft in your input file.'
   MSG_ERROR(message)
 end if



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
!O
!P 
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'scup_printbands',tread,'INT')
 if(tread==1) scup_dtset%scup_printbands=intarr(1)
 if(scup_dtset%scup_printbands<0 .or. scup_dtset%scup_printbands>1 )then
   write(message, '(a,I3,a,a,a,a,a)' )&
&   'scup_printbands is',scup_dtset%scup_printbands,', but the only allowed values',ch10,&
&   'are 0 and 1.',ch10,&
&   'Action: correct scup_printbands in your input file.'
   MSG_ERROR(message)
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'scup_printeigv',tread,'INT')
 if(tread==1) scup_dtset%scup_printeigv=intarr(1)
 if(scup_dtset%scup_printeigv<0 .or. scup_dtset%scup_printeigv>1 )then
   write(message, '(a,I3,a,a,a,a,a)' )&
&   'scup_printeigv is',scup_dtset%scup_printeigv,', but the only allowed values',ch10,&
&   'are 0 and 1.',ch10,&
&   'Action: correct scup_printeigv in your input file.'
   MSG_ERROR(message)
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'scup_printeltic',tread,'INT')
 if(tread==1) scup_dtset%scup_printeltic=intarr(1)
 if(scup_dtset%scup_printeltic<0 .or. scup_dtset%scup_printeltic>1 )then
   write(message, '(a,I3,a,a,a,a,a)' )&
&   'scup_printeltic is',scup_dtset%scup_printeltic,', but the only allowed values',ch10,&
&   'are 0 and 1.',ch10,&
&   'Action: correct scup_printeltic in your input file.'
   MSG_ERROR(message)
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'scup_printgeom',tread,'INT')
 if(tread==1) scup_dtset%scup_printgeom=intarr(1)
 if(scup_dtset%scup_printgeom<0 .or. scup_dtset%scup_printgeom>1 )then
   write(message, '(a,I3,a,a,a,a,a)' )&
&   'scup_printgeom is',scup_dtset%scup_printgeom,', but the only allowed values',ch10,&
&   'are 0 and 1.',ch10,&
&   'Action: correct scup_printgeom in your input file.'
   MSG_ERROR(message)
 end if

 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'scup_printorbocc',tread,'INT')
 if(tread==1) scup_dtset%scup_printorbocc=intarr(1)
 if(scup_dtset%scup_printorbocc<0 .or. scup_dtset%scup_printorbocc>1 )then
   write(message, '(a,I3,a,a,a,a,a)' )&
&   'scup_printorbocc is',scup_dtset%scup_printorbocc,', but the only allowed values',ch10,&
&   'are 0 and 1.',ch10,&
&   'Action: correct scup_printorbocc in your input file.'
   MSG_ERROR(message)
 end if


!Q
!R
!S 
!T 
!U 
!V 
!W
!X
!Y
!Z 

end subroutine invars10scup


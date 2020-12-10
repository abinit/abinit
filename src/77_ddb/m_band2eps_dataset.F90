!!****m* ABINIT/m_band2eps_dataset
!! NAME
!!  m_band2eps_dataset
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

module m_band2eps_dataset

 use defs_basis
 use m_abicore
 use m_errors

 use m_parser,    only : intagm

 implicit none

 private

 public :: band2eps_dataset_type
 public :: band2eps_dtset_free
 public :: outvars_band2eps
 public :: invars11
!!***

!----------------------------------------------------------------------

!!****t* m_band2eps_dataset/band2eps_dataset_type
!! NAME
!! band2eps_dataset_type
!!
!! FUNCTION
!! The band2eps_dataset_type structured datatype
!! gather all the input variables for the band2eps code.
!!
!! SOURCE

 type band2eps_dataset_type

! Variables should be declared on separated lines in order to reduce the occurence of git conflicts.
! Since all these input variables are described in the band2eps_help.html
! file, they are not described in length here ...
! Integer
  integer :: natom
  integer :: cunit
  integer :: nlines
  integer :: ngrad
  integer :: prtout
  real(dp) :: min
  real(dp) :: max


  integer,allocatable :: nqline(:)
! nqline(nlines)

  real(dp),allocatable :: scale(:)
! scale(nlines)

  integer,allocatable :: red(:)
! red(natom)

  integer,allocatable :: blue(:)
! blue(natom)

  integer,allocatable :: green(:)
! green(natom)

  character(len=6),allocatable :: qpoint_name(:)
! qpoint_name(nlines+1)


 end type band2eps_dataset_type
!!***

contains
!!***

!!****f* m_band2eps_dataset/band2eps_dtset_free
!!
!! NAME
!!   band2eps_dtset_free
!!
!! FUNCTION
!!   deallocate remaining arrays in the band2eps_dtset datastructure
!!
!! INPUTS
!!  band2eps_dtset = band2eps datastructure
!!
!! PARENTS
!!      band2eps
!!
!! CHILDREN
!!
!! NOTES
!!
!! SOURCE

subroutine band2eps_dtset_free(band2eps_dtset)

!Arguments ------------------------------------
!scalars
 type(band2eps_dataset_type), intent(inout) :: band2eps_dtset

! *************************************************************************

 ABI_SFREE(band2eps_dtset%nqline)
 ABI_SFREE(band2eps_dtset%scale)
 ABI_SFREE(band2eps_dtset%red)
 ABI_SFREE(band2eps_dtset%blue)
 ABI_SFREE(band2eps_dtset%green)
 ABI_SFREE(band2eps_dtset%qpoint_name)

end subroutine band2eps_dtset_free
!!***

!----------------------------------------------------------------------

!!****f* m_band2eps_dataset/invars11
!!
!! NAME
!! invars11
!!
!! FUNCTION
!! Open input file for the band2eps code, then reads or echoes the input information.
!!
!! INPUTS
!! lenstr=actual length of string
!! natom=number of atoms, needed for atifc
!! string*(*)=string of characters containing all input variables and data
!!
!! OUTPUT
!! band2eps_dtset= (derived datatype) contains all the input variables
!!
!! NOTES
!! Should be executed by one processor only.
!!
!!
!! PARENTS
!!      band2eps
!!
!! CHILDREN
!!
!! SOURCE

subroutine invars11 (band2eps_dtset,lenstr,string)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: lenstr
 character(len=strlen),intent(in) :: string
 type(band2eps_dataset_type),intent(inout) :: band2eps_dtset

!Local variables -------------------------
!Dummy arguments for subroutine 'intagm' to parse input file
!Set routine version number here:
!scalars
 integer :: ii,position,jdtset,marr,tread
 character(len=500) :: message
 character(len=fnlen) :: name_qpoint
!arrays
 integer,allocatable :: intarr(:)
 real(dp),allocatable :: dprarr(:)

!*********************************************************************
 marr=3
 ABI_ALLOCATE(intarr,(marr))
 ABI_ALLOCATE(dprarr,(marr))

 jdtset=1

! Need to get the number of atom
 band2eps_dtset%natom = 0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'natom',tread,'INT')
 if(tread==1) band2eps_dtset%natom=intarr(1)
 if (band2eps_dtset%natom <= tol10) then
   write(message,'(a,I5,3a)' )&
&   'natom ',band2eps_dtset%natom,', is not allowed ',ch10,&
&   'Action: correct natom in your input file.'
   ABI_ERROR(message)
 end if

! Need to get the number of lines
 band2eps_dtset%nlines = 0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'nlines',tread,'INT')
 if(tread==1) band2eps_dtset%nlines=intarr(1)
 if (band2eps_dtset%nlines <= tol10) then
   write(message,'(a,I5,3a)' )&
&   'nlines ',band2eps_dtset%nlines,', is not allowed ',ch10,&
&   'Action: correct nlines in your input file.'
   ABI_ERROR(message)
 end if

! Need to get the number of lines
 band2eps_dtset%cunit = 0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'cunit',tread,'INT')
 if(tread==1) band2eps_dtset%cunit=intarr(1)
 if (band2eps_dtset%cunit < 1 .or. band2eps_dtset%cunit > 2) then
   write(message,'(a,I5,3a)' )&
&   'cunit ',band2eps_dtset%cunit,', is not allowed ',ch10,&
&   'Action: correct cunit in your input file (1 or 2).'
   ABI_ERROR(message)
 end if

! Need to get the number of lines
 band2eps_dtset%ngrad = 0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'ngrad',tread,'INT')
 if(tread==1) band2eps_dtset%ngrad=intarr(1)
 if (band2eps_dtset%ngrad < 0) then
   write(message,'(a,I5,3a)' )&
&   'ngrad ',band2eps_dtset%ngrad,', is not allowed ',ch10,&
&   'Action: correct ngrad in your input file (positive value).'
   ABI_ERROR(message)
 end if

! Need to get the number of lines
 band2eps_dtset%min = zero
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'min',tread,'DPR')
 if(tread==1) band2eps_dtset%min=dprarr(1)

! Need to get the number of lines
 band2eps_dtset%max = zero
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'max',tread,'DPR')
 if(tread==1) band2eps_dtset%max=dprarr(1)

 ABI_ALLOCATE(band2eps_dtset%red,(band2eps_dtset%natom))
 ABI_ALLOCATE(band2eps_dtset%blue,(band2eps_dtset%natom))
 ABI_ALLOCATE(band2eps_dtset%green,(band2eps_dtset%natom))
 band2eps_dtset%red(:) = 0
 band2eps_dtset%blue(:) = 0
 band2eps_dtset%green(:) = 0

! Read prtout
 band2eps_dtset%prtout = 0
 call intagm(dprarr,intarr,jdtset,marr,1,string(1:lenstr),'prtout',tread,'INT')
 if(tread==1) band2eps_dtset%prtout=intarr(1)
 if (band2eps_dtset%prtout < 0 .or. band2eps_dtset%prtout > 1) then
   write(message,'(a,I5,3a)' )&
&   'prtout ',band2eps_dtset%prtout,', is not allowed ',ch10,&
&   'Action: correct prtout in your input file (0 or 1).'
   ABI_ERROR(message)
 end if

!natom dimension
 if(band2eps_dtset%natom > marr)then
   marr = band2eps_dtset%natom
   ABI_DEALLOCATE(intarr)
   ABI_DEALLOCATE(dprarr)
   ABI_ALLOCATE(intarr,(marr))
   ABI_ALLOCATE(dprarr,(marr))
 end if
 call intagm(dprarr,intarr,jdtset,marr,band2eps_dtset%natom,string(1:lenstr),'red',tread,'INT')
 if(tread==1) band2eps_dtset%red(1:band2eps_dtset%natom) = intarr(1:band2eps_dtset%natom)

 call intagm(dprarr,intarr,jdtset,marr,band2eps_dtset%natom,string(1:lenstr),'blue',tread,'INT')
 if(tread==1) band2eps_dtset%blue(1:band2eps_dtset%natom) = intarr(1:band2eps_dtset%natom)

 call intagm(dprarr,intarr,jdtset,marr,band2eps_dtset%natom,string(1:lenstr),'green',tread,'INT')
 if(tread==1) band2eps_dtset%green(1:band2eps_dtset%natom) = intarr(1:band2eps_dtset%natom)

!nlines dimenstion
 ABI_ALLOCATE(band2eps_dtset%nqline,(band2eps_dtset%nlines))
 ABI_ALLOCATE(band2eps_dtset%scale,(band2eps_dtset%nlines))
 band2eps_dtset%nqline(:) = 0
 band2eps_dtset%scale(:) = zero

 if(band2eps_dtset%nlines > marr)then
   marr = band2eps_dtset%nlines
   ABI_DEALLOCATE(intarr)
   ABI_DEALLOCATE(dprarr)
   ABI_ALLOCATE(intarr,(marr))
   ABI_ALLOCATE(dprarr,(marr))
 end if

 call intagm(dprarr,intarr,jdtset,marr,band2eps_dtset%nlines,string(1:lenstr),'nqline',tread,'INT')
 if(tread==1) band2eps_dtset%nqline(1:band2eps_dtset%nlines) = intarr(1:band2eps_dtset%nlines)
!REMOVE THE LAST LINE OF NQLINE
  band2eps_dtset%nqline(band2eps_dtset%nlines) = band2eps_dtset%nqline(band2eps_dtset%nlines) - 1

 call intagm(dprarr,intarr,jdtset,marr,band2eps_dtset%nlines,string(1:lenstr),'scale',tread,'DPR')
 if(tread==1) band2eps_dtset%scale(1:band2eps_dtset%nlines) = dprarr(1:band2eps_dtset%nlines)

!nline+1 dimension
 ABI_ALLOCATE(band2eps_dtset%qpoint_name,(band2eps_dtset%nlines+1))
 band2eps_dtset%qpoint_name(:) = ""
 position = index(string(1:lenstr),trim("QPOINT_NAME")) + 11
 name_qpoint = trim(string(position:(position + &
&          len(band2eps_dtset%qpoint_name(1))*band2eps_dtset%nlines+1)))
 read(name_qpoint,*) (band2eps_dtset%qpoint_name(ii),ii=1,band2eps_dtset%nlines+1)

end subroutine invars11
!!***

!----------------------------------------------------------------------

!!****f* m_band2eps_dataset/outvars_band2eps
!!
!! NAME
!! outvars_band2eps
!!
!! FUNCTION
!! Open input file for the band2eps code, then
!! echoes the input information.
!!
!! INPUTS
!! band2eps_dtset= (derived datatype) contains all the input variables
!! nunit=unit number for input or output
!!
!! OUTPUT
!!  (only writing)
!!
!! NOTES
!! Should be executed by one processor only.
!!
!! PARENTS
!!      band2eps
!!
!! CHILDREN
!!
!! SOURCE

subroutine outvars_band2eps (band2eps_dtset,nunit)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: nunit
 type(band2eps_dataset_type),intent(in) :: band2eps_dtset

!Local variables -------------------------
!Set routine version number here:
!scalars
 integer :: ii

!*********************************************************************

!Write the heading
 write(nunit,'(80a,a)') ('=',ii=1,80),ch10
 write(nunit, '(a,a)' )&
& ' -outvars_band2eps: echo values of input variables ----------------------',ch10

!The flags
 if(band2eps_dtset%natom/=0 )then
   write(nunit,'(a)')' Informations :'
   if(band2eps_dtset%natom/=0)write(nunit,'(3x,a9,3i10)')'  natom',band2eps_dtset%natom
   if(band2eps_dtset%cunit/=0)write(nunit,'(3x,a9,3i10)')'  cunit',band2eps_dtset%cunit
   if(band2eps_dtset%nlines/=0)write(nunit,'(3x,a9,3i10)')' nlines',band2eps_dtset%nlines
   if(band2eps_dtset%min/=0)write(nunit,'(3x,a9,3f12.4)')'    Min',band2eps_dtset%min
   if(band2eps_dtset%max/=0)write(nunit,'(3x,a9,3f12.4)')'    Max',band2eps_dtset%max
   if(band2eps_dtset%ngrad/=0)write(nunit,'(3x,a9,3i10)')'  ngrad',band2eps_dtset%ngrad
   write(nunit,'(3x,a9,6i10)')'    red',(band2eps_dtset%red(ii),ii=1,band2eps_dtset%natom)
   write(nunit,'(3x,a9,6i10)')'   blue',(band2eps_dtset%blue(ii),ii=1,band2eps_dtset%natom)
   write(nunit,'(3x,a9,6i10)')'  green',(band2eps_dtset%green(ii),ii=1,band2eps_dtset%natom)
   write(nunit,'(3x,a9,8i10)')' nqline',(band2eps_dtset%nqline(ii),ii=1,band2eps_dtset%nlines)
   write(nunit,'(3x,a9,8a)')'    point',(band2eps_dtset%qpoint_name(ii),ii=1,band2eps_dtset%nlines+1)
 end if

 write(nunit,'(a,80a,a)') ch10,('=',ii=1,80),ch10

end subroutine outvars_band2eps
!!***

!----------------------------------------------------------------------

end module m_band2eps_dataset
!!***

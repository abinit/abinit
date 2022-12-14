  !!****m* ABINIT/m_multibinit_main2
  !! NAME
  !! m_multibinit_main2
  !!
  !! FUNCTION
  !! multibinit_main2 is the main function of multibinit. It runs after
  !! the filenames from files file is read. It then does everything else.
  !!
  !! Datatypes:
  !!
  !!
  !! Subroutines:
  !! multibinit_main2
  !!
  !!
  !! COPYRIGHT
  !! Copyright (C) 2001-2022 ABINIT group (hexu)
  !! This file is distributed under the terms of the
  !! GNU General Public License, see ~abinit/COPYING
  !! or http://www.gnu.org/copyleft/gpl.txt .
  !! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
  !!
  !! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"
module m_multibinit_main2
  use defs_basis
  use defs_abitypes
  use m_build_info
  use m_xmpi
  use m_xomp
  use m_abicore
  use m_errors
  use m_multibinit_manager, only: mb_manager_t
  
  implicit none

  !!***
contains

  !!****f* m_multibinit_main2/multibinit_main2
  !!
  !! NAME
  !! multibinit_main2
  !!
  !! FUNCTION
  !! The main function
  !!
  !! INPUTS
  !! filnam: file names from files file
  !!
  !! OUTPUT
  !!
  !! SOURCE
  subroutine multibinit_main2(input_path, filnam, dry_run)
    character(len=fnlen), intent(inout) :: input_path
    character(len=fnlen), intent(inout) :: filnam(18)
    integer, intent(in) :: dry_run
    type(mb_manager_t) :: manager
    call manager%run_all(input_path, filnam, dry_run)
  end subroutine multibinit_main2
  !!***


!!****f* m_multibinit_main/herald_multibinit
!! NAME
!!  herald_multibinit
!!
!! FUNCTION
!!  Prints out a message to unit iout giving info about current
!!  code, version of code, platform, and starting date.
!!  Modified from m_special_msg.F90/herald
!!
!! INPUTS
!!  code_name= code name
!!  code_version= code version
!!  iout=unit number for output
!!
!! OUTPUT
!!  (only writing)
!!
!! SOURCE

subroutine herald_multibinit(code_name,code_version,iout, lattmode)

!Arguments ------------------------------------
 integer,intent(in) :: iout
 character(len=*),intent(in) :: code_name
 character(len=*),intent(in) :: code_version
 logical, intent(in) :: lattmode ! 1 lattice, 0 other

!Local variables-------------------------------
 integer :: day,dd,ja,jy,jm,jdn,mm,mm_rel,year,year_rel
 integer :: values(8)
 character(len=5) :: strzone
 character(len=8) :: strdat
 character(len=10) :: strtime
 character(len=500) :: msg
 character(len=3),parameter :: day_names(7)=(/'Mon','Tue','Wed','Thu','Fri','Sat','Sun'/)
 character(len=3),parameter :: month_names(12)=(/'Jan','Feb','Mar','Apr','May','Jun',&
                                                 'Jul','Aug','Sep','Oct','Nov','Dec'/)

! *************************************************************************

!RELEASE TIME FROM ABIRULES
 year_rel=2022
 mm_rel=07
!END OF RELEASE TIME

 write(iout, '(a,/,a,/,a,/,a)' ) &
'******************************************************************************************', &
'                                Welcome to MULTIBINIT,                         ', &
' a software platform designed for the construction and use of second-principles models', &
'                   for lattice, spin and electron degrees of freedom.'


!The technique used hereafter is the only one that we have found to obtain
!perfect transferability across platforms and OS.
 write(iout, '(/,a,a,a,a,a)' ) '.Version ',trim(code_version),' of ',trim(code_name),' '
#if defined HAVE_MPI
 write(iout, '(a,a,a,/)' ) '.(MPI version, prepared for a ',build_target,' computer) '
#else
 write(iout, '(a,a,a,/)' ) '.(sequential version, prepared for a ',build_target,' computer) '
#endif



!GNU GPL license
 write(iout, '(a,/,a,a,a,/,a,/,a,/,a,/)' ) &
 '.Copyright (C) 1998-2022 ABINIT group . ',&
 ' ',trim(code_name),' comes with ABSOLUTELY NO WARRANTY.',&
 ' It is free software, and you are welcome to redistribute it',&
 ' under certain conditions (GNU General Public License,',&
 ' see ~abinit/COPYING or http://www.gnu.org/copyleft/gpl.txt).'


if(.not. lattmode) then
 write(iout, '(a,/,a,/,a,/,a,/,a,/)' ) &
 ' ABINIT is a project of the Universite Catholique de Louvain,',&
 ' Corning Inc. and other collaborators, see ~abinit/doc/developers/contributors.txt .',&
 ' Please read https://docs.abinit.org/theory/acknowledgments for suggested',&
 ' acknowledgments of the ABINIT effort.',&
 ' For more information, see https://www.abinit.org .'
endif

 write(iout, '(a,/,a,/)' ) &
' MULTIBINIT is a software project of the University of Liege', &
' (PHYTHEMA & NANOMAT groups), in collaboration with other partners.'

 if (lattmode) then
   write(iout, '(a,/,/,a,/,/,a,/, a,/,/, a,/, a,/, a,/,/, a,/,/, a/,/)' ) &
'-----------------------------------------------------------------------------------------', &
'                          MULTIBINIT – LATTICE MODELS                   ', &
' Project initiated and coordinated by Philippe GHOSEZ and his group at ULiege', &
'   (Philippe.Ghosez@uliege.be).',  &
' Main contributors: Alexandre MARTIN, Jordan BIEDER, Michael Marcus SCHMITT,', &
'   Louis BASTOGNE, Xu HE, Alireza SASANI, Huazhang ZHANG, Subhadeep BANDYOPADHYAY,', &
'   Philippe GHOSEZ.', &
' Technical support: Xu HE (X.He@uliege.be)', &
'*****************************************************************************************'
 else
   write(iout, '(a,/)' ) &
   '*****************************************************************************************'
 end if



!Get year, month and day
 call date_and_time(strdat,strtime,strzone,values)
 year=values(1)
 mm=values(2)
 dd=values(3)

!Get day of the week
 if (mm.gt.2) then
   jy=year
   jm=mm+1
 else
   jy=year-1
   jm=mm+13
 end if
 jdn=int(365.25d0*jy)+int(30.6001d0*jm)+dd+1720995
 ja=int(0.01d0*jy)
 jdn=jdn+2-ja+int(quarter*ja)
 day=mod(jdn,7)+1

!Print date in nice format (* new format *)
 write(iout, '(/,a,a,1x,i2,1x,a,1x,i4,a,/,a,i2.2,a,i2.2,a)' ) &
 '.Starting date : ',day_names(day),dd,month_names(mm),year,'.','- ( at ',values(5),'h',values(6),' )'
 write(iout,*)' '

!Impose a maximal life cycle of 3 years
 if(year>year_rel+3 .or. (year==year_rel+3 .and. mm>mm_rel) ) then
   write(msg, '(5a,i4,5a)' )&
   '- The starting date is more than 3 years after the initial release',ch10,&
   '- of this version of ABINIT, namely ',month_names(mm_rel),' ',year_rel,'.',ch10,&
   '- This version of ABINIT is not supported anymore.',ch10,&
   '- Action: please, switch to a more recent version of ABINIT.'
   call wrtout(iout,msg,'COLL')

!  Gives a warning beyond 2 years
 else if(year>year_rel+2 .or. (year==year_rel+2 .and. mm>mm_rel) ) then
   write(msg, '(5a,i4,6a)' )&
   '- The starting date is more than 2 years after the initial release',ch10,&
   '- of this version of ABINIT, namely ',month_names(mm_rel),' ',year_rel,'.',ch10,&
   '- Note that the use beyond 3 years after the release will not be supported.',ch10,&
   '- Action: please, switch to a more recent version of ABINIT.',ch10
   call wrtout(iout,msg,'COLL')
 end if

end subroutine herald_multibinit
!!***



end module m_multibinit_main2

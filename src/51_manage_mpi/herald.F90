!{\src2tex{textfont=tt}}
!!****f* ABINIT/herald
!! NAME
!!  herald
!!
!! FUNCTION
!!  Prints out a message to unit iout giving info about current
!!  code, version of code, platform, and starting date.
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2016 ABINIT group (DCA, XG, GMR, LSI, MM, MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  code_name= code name
!!  code_version= code version
!!  iout=unit number for output
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      abinit,aim,anaddb,bsepostproc,cut3d,fftprof,ioprof,lapackprof,mrgddb
!!      mrgdv,mrggkk,mrgscr,optic,ujdet,vdw_kernelgen
!!
!! CHILDREN
!!      date_and_time,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine herald(code_name,code_version,iout)

 use defs_basis
 use m_build_info
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'herald'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: iout
 character(len=*),intent(in) :: code_name
 character(len=*),intent(in) :: code_version

!Local variables-------------------------------
 integer :: day,dd,ja,jy,jm,jdn,mm,mm_rel,year,year_rel
 integer :: values(8)
 character(len=5) :: strzone
 character(len=8) :: strdat
 character(len=10) :: strtime
 character(len=500) :: message
 character(len=3),parameter :: day_names(7)=(/'Mon','Tue','Wed','Thu','Fri','Sat','Sun'/)
 character(len=3),parameter :: month_names(12)=(/'Jan','Feb','Mar','Apr','May','Jun',&
&                                                'Jul','Aug','Sep','Oct','Nov','Dec'/)

! *************************************************************************

!RELEASE TIME FROM ABIRULES
 year_rel=2016
 mm_rel=11
!END OF RELEASE TIME

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
& '.Copyright (C) 1998-2016 ABINIT group . ',&
& ' ',trim(code_name),' comes with ABSOLUTELY NO WARRANTY.',&
& ' It is free software, and you are welcome to redistribute it',&
& ' under certain conditions (GNU General Public License,',&
& ' see ~abinit/COPYING or http://www.gnu.org/copyleft/gpl.txt).'

 if(trim(code_name)=='OPTIC')then
   write(iout, '(a,a,a,/,a,/,a,/,a,/,a,/)' ) &
&   ' ',trim(code_name),' has originally been developed by',&
&   ' Sangeeta Sharma and incorporated in ABINIT with the help of M. Verstraete.',&
&   ' Please refer to : ',&
&   ' S. Sharma, J. K. Dewhurst and C. Ambrosch-Draxl, Phys. Rev. B 67, 165332 (2003), and',&
&   ' S. Sharma and C. Ambrosch-Draxl, Physica Scripta T 109 (2004).'
 end if

 write(iout, '(a,/,a,/,a,/,a,/,a)' ) &
& ' ABINIT is a project of the Universite Catholique de Louvain,',&
& ' Corning Inc. and other collaborators, see ~abinit/doc/developers/contributors.txt .',&
& ' Please read ~abinit/doc/users/acknowledgments.html for suggested',&
& ' acknowledgments of the ABINIT effort.',&
& ' For more information, see http://www.abinit.org .'

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
& '.Starting date : ',day_names(day),dd,month_names(mm),year,'.','- ( at ',values(5),'h',values(6),' )'
 write(iout,*)' '

!Impose a maximal life cycle of 3 years
 if(year>year_rel+3 .or. (year==year_rel+3 .and. mm>mm_rel) ) then 
   write(message, '(5a,i4,5a)' )&
&   '- The starting date is more than 3 years after the initial release',ch10,&
&   '- of this version of ABINIT, namely ',month_names(mm_rel),' ',year_rel,'.',ch10,&
&   '- This version of ABINIT is not supported anymore.',ch10,&
&   '- Action: please, switch to a more recent version of ABINIT.'
   call wrtout(iout,message,'COLL')

!  Gives a warning beyond 2 years
 else if(year>year_rel+2 .or. (year==year_rel+2 .and. mm>mm_rel) ) then
   write(message, '(5a,i4,6a)' )&
&   '- The starting date is more than 2 years after the initial release',ch10,&
&   '- of this version of ABINIT, namely ',month_names(mm_rel),' ',year_rel,'.',ch10,&
&   '- Note that the use beyond 3 years after the release will not be supported.',ch10,&
&   '- Action: please, switch to a more recent version of ABINIT.',ch10
   call wrtout(iout,message,'COLL')
 end if

end subroutine herald
!!***

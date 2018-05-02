!{\src2tex{textfont=tt}}
!!****f* ABINIT/status
!! NAME
!! status
!!
!! FUNCTION
!! Routine for description of the status of the calculation
!! Eventually open the status file, write different information,
!! and close the file. The output rate and shift are governed by istat
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (XG,TD)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  counter=value of the loop counter at that level
!!  file (optional argument)=name of the status file
!!  istat=gives the rate or shift of output. The status file will be opened
!!     and written only once every "istatr" calls.
!!     This variable is saved at the first call (just after the first
!!     call to invars0. The shift "istatshft" is saved at the second call.
!!     In subsequent calls, istat has no meaning.
!!  level=number of the level of the calling subroutine (see the description later)
!!  routine=string of 14 characters indicating the status inside the level
!!
!! OUTPUT
!!  (only writing)
!!
!! NOTES
!! Warning : The string "routine" can have any size but
!! it is truncated to a size of 14.
!! because of the behaviour of some compilers, the string
!! "routine" should always have 14 characters in the calling subroutine
!!
!! PARENTS
!!      abinit,dfpt_accrho,dfpt_looppert,dfpt_rhofermi,dfpt_scfcv,dfpt_vtowfk
!!      dfpt_wfkfermi,dfptnl_loop,dfptnl_mv,dfptnl_resp,driver,gstate,gstateimg
!!      m_ab7_invars_f90,mover,nonlinear,respfn,scfcv,vtorhotf
!!
!! CHILDREN
!!      timab
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine status(counter,filstat,istat,level,routine)

 use defs_basis
 use m_errors

 use m_io_tools,   only : open_file
 use m_time,       only : timab

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'status'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: counter,istat,level
 character(len=*),intent(in) :: routine
 character(len=*),intent(in) :: filstat

!Local variables-------------------------------
!scalars
 integer,parameter :: mcounter=2,mlevel=120
 integer,save :: output_rate=1,shift_rate=1,statnu=0
 integer :: ilevel,temp_unit
 character(len=12) :: headwr
 character(len=500) :: message
!arrays
 integer,save :: active(mlevel),actual_counter(mlevel,mcounter)
 integer,save :: ncounter(mlevel)
 integer,save :: list_level(29)=&
&  (/1,2,3,100,101,102,110,111,112,113,114,10,11,12,13,14,15,16,17,18,20,30,31,40,41,50,51,52,90/)
 real(dp) :: tsec(2)
 character(len=20),save :: nm_levels(mlevel),nm_routine(mlevel)
 character(len=12),save :: nm_counter(mlevel,mcounter)

!***********************************************************************

 if (.not.do_write_status .or. output_rate==0) return

 call timab(73,1,tsec)

!Note : all processors have their own file, so no special
!attention must be paid to the parallel case.
!Initialisation
 if(statnu==0)then
   nm_routine(:)='                  '
   active(:)=0
   actual_counter(:,:)=0

!  List of names for each level
!  Numbers from 1 to 9 are for abinit and driver
!  Numbers from 100 to 120 are for optdriver=0 routines (GS)
!  Numbers between 10 and 19 are for optdriver=1 routines (RF)
!  Numbers between 20 and 29 are for optdriver=2 routines (suscep)
!  Numbers between 30 and 39 are for optdriver=3 routines (screening)
!  Numbers between 40 and 49 are for optdriver=4 routines (sigma)
!  Numbers between 50 and 59 are for optdriver=5 routines (nonlinear)
!  When you add a level number, or modify one, do not forget to change list_level

   nm_levels(1)   ='abinit              '
   ncounter(1)=0
   nm_counter(1,1)='            '

   nm_levels(2)   ='driver              '
   ncounter(2)=1
   nm_counter(2,1)='jdtset     ='

   nm_levels(3)   ='ab7_invars_load     '
   ncounter(3)=0
   nm_counter(3,1)='            '

!  Optdriver=0
   nm_levels(100)   ='gstateimg           '
   ncounter(100)=2
   nm_counter(100,1)='idynimage  ='
   nm_counter(100,2)='itimimage  ='

   nm_levels(101)   ='gstate              '
   ncounter(101)=1
   nm_counter(101,1)='itime      ='

   nm_levels(102)   ='mover               '
   ncounter(102)=2
   nm_counter(102,1)='icycle     ='
   nm_counter(102,2)='itime      ='

   nm_levels(110)   ='scfcv               '
   ncounter(110)=1
   nm_counter(110,1)='istep      ='

   nm_levels(111)   ='vtorho(tf)          '
   ncounter(111)=2
   nm_counter(111,1)='isppol     ='
   nm_counter(111,2)='ikpt       ='

   nm_levels(112)   ='vtowfk              '
   ncounter(112)=2
   nm_counter(112,1)='inonsc     ='
   nm_counter(112,2)='iband      ='

   nm_levels(113)   ='cgwf                '
   ncounter(113)=1
   nm_counter(113,1)='iline      ='

   nm_levels(114)   ='getghc              '
   ncounter(114)=0
   nm_counter(114,1)='            '


!  Optdriver=1
   nm_levels(10)   ='respfn              '
   ncounter(10)=0
   nm_counter(10,1)='            '

   nm_levels(11)   ='dfpt_looppert       '
   ncounter(11)=1
   nm_counter(11,1)='respcase   ='

   nm_levels(12)   ='dfpt_scfcv          '
   ncounter(12)=1
   nm_counter(12,1)='istep      ='

   nm_levels(13)   ='dfpt_vtorho         '
   ncounter(13)=2
   nm_counter(13,1)='isppol     ='
   nm_counter(13,2)='ikpt       ='

   nm_levels(14)   ='dfpt_vtowfk         '
   ncounter(14)=2
   nm_counter(14,1)='inonsc     ='
   nm_counter(14,2)='iband      ='

   nm_levels(15)   ='dfpt_cgwf           '
   ncounter(15)=1
   nm_counter(15,1)='iline      ='

   nm_levels(16)   ='getgh1c             '
   ncounter(16)=0
   nm_counter(16,1)='            '

   nm_levels(17)   ='dfpt_rhofermi       '
   ncounter(17)=2
   nm_counter(17,1)='isppol     ='
   nm_counter(17,2)='ikpt       ='

   nm_levels(18)   ='dfpt_wfkfermi       '
   ncounter(18)=2
   nm_counter(18,1)='inonsc     ='
   nm_counter(18,2)='iband      ='


!  Optdriver=2
   nm_levels(20)   ='suscep              '
   ncounter(20)=0
   nm_counter(20,1)='            '


!  Optdriver=3
   nm_levels(30)   ='screening           '
   ncounter(30)=1
   nm_counter(30,1)='iqpt       ='

!  Optdriver=4
   nm_levels(40)   ='sigma               '
   ncounter(40)=1
   nm_counter(40,1)='ikpt_gw    ='

!  Optdriver=5
   nm_levels(50)   ='nonlinear           '
   ncounter(50)=0
   nm_counter(50,1)='            '

   nm_levels(51)   ='dfptnl_loop         '
   ncounter(51)=2
   nm_counter(51,1)='pert1case  ='
   nm_counter(51,2)='pert3case  ='

   nm_levels(52)   ='mv_/dfptnl_resp     '
   ncounter(52)=2
   nm_counter(52,2)='ikpt       ='

!  Optdriver=9
   nm_levels(90)   ='bethe_salpether     '
   ncounter(90)=0
   nm_counter(90,1)='            '

   if(istat<0)then
     write(message, '(a,i7,a,a,a,a,a)' )&
&     'the value of the input variable istatr is',istat,' ,',ch10,&
&     'while it must be a positive number.',ch10,&
&     'Action : change istatr in your input file.'
     MSG_ERROR(message)
   end if
   output_rate=istat
 end if

!The true input variable "shift_rate" is only available at the second call
 if(statnu==1)then
   if(istat<0 .or. istat>output_rate)then
     write(message, '(a,i7,3a,i7,2a)' )&
&     'the value of the input variable istatshft is',istat,' ,',ch10,&
&     'while it must be a positive number smaller or equal to istatr=',output_rate,ch10,&
&     'Action: change istatshft in your input file.'
     MSG_ERROR(message)
   end if
   shift_rate=istat
   if(shift_rate==output_rate)shift_rate=0

!  At this second call, also feed information that the abinit routine called ab7_invars_load
   write(unit=nm_routine(1),fmt='(a20)') 'call ab7_invars_load'
   active(1)=1
 end if

!Check the value of level
 if( minval(abs(list_level(:)-level)) /= 0)then
   write(message, '(a,i5,3a)' )&
&   '  The value of level in the calling routine is',level,' ,',ch10,&
&   '  which is not an allowed value.'
   MSG_BUG(message)
 end if

!Assign the info about the actual routine
 write(unit=nm_routine(level),fmt='(a20)') trim(adjustl(routine))
 if(trim(adjustl(nm_routine(level)))=='exit')then
!  The value 2 will be changed to 0 at the end of the routine.
   active(level)=2
 else if(trim(adjustl(nm_routine(level)))=='')then
   active(level)=0
 else
   active(level)=1
 end if

!Assign the info about the actual counter
 if(counter>=0)then
   if(ncounter(level)==1)then
     actual_counter(level,1)=counter
   else if(ncounter(level)==2)then
     actual_counter(level,2)=counter/100
!    The counter number 1 does not allow more than 99 passes
     actual_counter(level,1)=counter-(counter/100)*100
   end if
 end if

!============================================================

!After treatment of present information, output of the status
 statnu=statnu+1

!DEBUG
! write(std_out,*)' status : statnu, output_rate, shift_rate=',statnu,output_rate, shift_rate
! write(std_out,*)'level,routine=',level,routine
! write(std_out,*)'active(level)=',active(level)
! write(std_out,*)'counter,actual_counter(level,1:2)=',counter,actual_counter(level,1:2)
! write(std_out,*)'List of active levels :'
! do ilevel=1,mlevel
!   if(active(ilevel)/=0)write(std_out,*)' Active level number=',ilevel
! end do
!ENDDEBUG

 if(statnu > 2 )then
   if( mod(statnu,output_rate)==shift_rate )then

     if (open_file(filstat,message,newunit=temp_unit,form='formatted',status='unknown') /= 0) then
       MSG_ERROR(message)
     end if

     rewind temp_unit
     write(temp_unit,*)

     headwr='(a,i4,a,i6 )'
     if(statnu>=100000)   headwr='(a,i4,a,i9 )'
     if(output_rate>=1000)headwr='(a,i6,a,i6 )'
     if(statnu>=100000 .and. output_rate>=1000)   headwr='(a,i6,a,i9 )'
     if(statnu>=100000000)headwr='(a,i6,a,i12)'
     write(temp_unit,headwr)' Status file, with repetition rate',output_rate,', status number',statnu
     write(temp_unit,*)

!    Treat every level one after the other
     do ilevel=1,mlevel
!      This level must be activated in order to have a corresponding output
       if(active(ilevel)==1 .or. active(ilevel)==2)then

         write(temp_unit,'(4a)')'  Level ',nm_levels(ilevel),' : ',adjustl(nm_routine(ilevel))

!        Is there a counter for this level ?
         if(ncounter(ilevel)>=1)then

           if(actual_counter(ilevel,1)>0)then
             write(temp_unit,'(a,a,i5)')'  ',nm_counter(ilevel,1),actual_counter(ilevel,1)
           end if
           if(ncounter(ilevel)==2)then
             if(actual_counter(ilevel,2)>0)then
               write(temp_unit,'(a,a,i5)')'  ',nm_counter(ilevel,2),actual_counter(ilevel,2)
             end if
           end if

         end if
       end if ! End of the check on activation of the level
     end do ! End of the loop on the levels

     close (temp_unit)
   end if ! End of the repetition rate check
 end if ! statnu > 2

 if (active(level)==2)then
   active(level)=0
   nm_routine(level)='              '
 end if

 call timab(73,2,tsec)

end subroutine status
!!***

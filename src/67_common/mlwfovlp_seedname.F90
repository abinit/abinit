!{\src2tex{textfont=tt}}
!!****f* ABINIT/mlwfovlp_seedname
!! NAME
!! mlwfovlp_seedname
!!
!! FUNCTION
!! Get seed name and file names of all wannier90 related files
!!
!! COPYRIGHT
!!  Copyright (C) 2009-2016 ABINIT group (T Rangel, DRH)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! fname_w90=root name of file appended with _w90
!!
!! OUTPUT
!! filew90_win= main input file for Wannier90
!! filew90_wout= main output file for Wannier90
!! filew90_amn= file containing Amn matrix
!! filew90_ramn= file containing Amn matrix (random initial projections)
!! filew90_mmn= file containing Mmn matrix
!! filew90_eig= file containing eigenvalues
!! nsppol= number of spin polarizations
!! seed_name= common seed name for all wannier90 related files
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      mlwfovlp
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine mlwfovlp_seedname(fname_w90,filew90_win,filew90_wout,filew90_amn,&
& filew90_ramn,filew90_mmn,filew90_eig,nsppol,seed_name)
    
 use defs_basis
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mlwfovlp_seedname'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: nsppol
 character(len=fnlen),intent(out) :: filew90_win(nsppol),filew90_wout(nsppol),filew90_amn(nsppol),filew90_ramn(nsppol)
 character(len=fnlen),intent(out) :: filew90_mmn(nsppol),filew90_eig(nsppol),seed_name(nsppol)
 character(len=fnlen),intent(in) :: fname_w90

!Local variables-------------------------------
 integer::isppol
 character(len=fnlen) :: test_win1,test_win2,test_win3
 logical :: lfile
 character(len=500) :: message                   ! to be uncommented, if needed
 character(len=10)::postfix 
! *************************************************************************
 
 seed_name(:)=trim(fname_w90)
 do isppol=1,nsppol
   if(nsppol==1)postfix='.win'
   if(nsppol==2 .and. isppol==1)postfix='_up.win'
   if(nsppol==2 .and. isppol==2)postfix='_down.win'
!  
   filew90_win(isppol)=trim(seed_name(isppol))//trim(postfix)
   test_win1=filew90_win(isppol)
   inquire(file=filew90_win(isppol),exist=lfile)

   if(.not.lfile) then
     seed_name(isppol)='wannier90'
     filew90_win(isppol)=trim(seed_name(isppol))//trim(postfix)
     test_win2=filew90_win(isppol)
     inquire(file=filew90_win(isppol),exist=lfile)
   end if

   if(.not.lfile) then
     seed_name(isppol)='w90'
     filew90_win=trim(seed_name(isppol))//trim(postfix)
     test_win3=filew90_win(isppol)
     inquire(file=filew90_win(isppol),exist=lfile)
   end if
   
   if(.not.lfile) then
     write(message,'(17a)')ch10,&
&     ' mlwfovlp_seedname : ERROR - ',ch10,&
&     ' wannier90 interface needs one of the following files:',ch10,&
&     '      ',trim(test_win1),ch10,&
&     '      ',trim(test_win2),ch10,&
&     '      ',trim(test_win3),ch10,&
&     ' Action: read wannier90 tutorial and/or user manual',ch10,&
&     '  and supply proper *.win file'
     MSG_ERROR(message)
   end if
 end do !isppol


!Files having different names for
!different spin polarizations
 if(nsppol==1) then
   filew90_win(1) =trim(seed_name(1))//'.win'
   filew90_wout(1)=trim(seed_name(1))//'.wout'
   filew90_ramn(1)=trim(seed_name(1))//'random.amn'
   filew90_amn(1) =trim(seed_name(1))//'.amn'
   filew90_mmn(1) =trim(seed_name(1))//'.mmn'
   filew90_eig(1) =trim(seed_name(1))//'.eig'
 elseif(nsppol==2) then
   filew90_win(1) =trim(seed_name(1))//'_up.win'
   filew90_win(2) =trim(seed_name(2))//'_down.win'
!  
   filew90_wout(1)=trim(seed_name(1))//'_up.wout'
   filew90_wout(2)=trim(seed_name(2))//'_down.wout'
!  
   filew90_ramn(1)=trim(seed_name(1))//'random_up.amn'
   filew90_ramn(2)=trim(seed_name(2))//'random_down.amn'
!  
   filew90_amn(1)=trim(seed_name(1))//'_up.amn'
   filew90_amn(2)=trim(seed_name(2))//'_down.amn'
!  
   filew90_mmn(1)=trim(seed_name(1))//'_up.mmn'
   filew90_mmn(2)=trim(seed_name(2))//'_down.mmn'
!  
   filew90_eig(1)=trim(seed_name(1))//'_up.eig'
   filew90_eig(2)=trim(seed_name(2))//'_down.eig'
 end if
!change also seed_name for nsppol=2
 if(nsppol==2) then
   seed_name(1)=trim(seed_name(1))//'_up'
   seed_name(2)=trim(seed_name(2))//'_down'
 end if
!End file-name section

 write(message, '(a,a)' ) ch10,&
& '---------------------------------------------------------------'
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')
 write(message, '(5a)' ) ch10,&
& '  Calculation of overlap and call to wannier90 library ',ch10,&
& '  to obtain maximally localized wannier functions ',ch10

 call wrtout(std_out,  message,'COLL')
 call wrtout(ab_out,message,'COLL')

 if(nsppol==1) then
   write(message, '(23a)' ) &
&   '  - ',trim(filew90_win(1)),' is a mandatory secondary input',ch10,&
&   '  - ',trim(filew90_wout(1)),' is the output for the library',ch10,&
&   '  - ',trim(filew90_ramn(1)),' contains random projections',ch10,&
&   '  - ',trim(filew90_amn(1)),' contains projections',ch10,&
&   '  - ',trim(filew90_mmn(1)),' contains the overlap',ch10,&
&   '  - ',trim(filew90_eig(1)),' contains the eigenvalues'
 elseif(nsppol==2) then
   write(message, '(41a)' ) & 
&   '  - ',trim(filew90_win(1)),&
&   ' and ',trim(filew90_win(2)),ch10,'are mandatory secondary input',ch10,&
&   '  - ',trim(filew90_wout(1)),&
&   ' and ',trim(filew90_wout(2)),ch10,' are the output for the library',ch10,&
&   '  - ',trim(filew90_ramn(1)),&
&   ' and ',trim(filew90_ramn(2)),ch10,' contain random projections',ch10,&
&   '  - ',trim(filew90_amn(1)),&
&   ' and ',trim(filew90_amn(2)),ch10,' contain projections',ch10,&
&   '  - ',trim(filew90_mmn(1)),&
&   ' and ',trim(filew90_mmn(2)),ch10,' contain the overlap',ch10,&
&   '  - ',trim(filew90_eig(1)),&
&   ' and ',trim(filew90_eig(2)),ch10,' contain the eigenvalues'
 end if
 call wrtout(std_out,  message,'COLL')
 call wrtout(ab_out,message,'COLL')

 write(message, '(a,a)' ) ch10,&
& '---------------------------------------------------------------'
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')

end subroutine mlwfovlp_seedname
!!***



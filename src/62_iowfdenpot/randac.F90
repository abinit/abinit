!{\src2tex{textfont=tt}}
!!****f* ABINIT/randac
!! NAME
!! randac
!!
!! FUNCTION
!! Random access to a wavefunction file : select the proper record in the wavefunction file.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! This routine originally written by Zachary Levine
!!
!! INPUTS
!!  debug=if > 0, prints debugging output
!!   without change the reading would begin at the wavefunction block ikptsp_prev+1 .
!!  headform1=format of the header of the wf file, also needed for reading the block
!!  ikpt=the k point at which data is desired
!!  isppol=spin channel
!!  nband(nkpt*nsppol)=number of bands at each k point and spin polarization
!!  nkpt=number of k points
!!  nsppol=number of spin channels
!!  wffinp=structured info for reading the wavefunction
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  ikptsp_prev=at input:
!!     number of the previously read k point - spin polarisation; 0 if new file;
!!     without change the reading would begin at the wavefunction block ikptsp_prev+1 .
!!    at output: value of ikptsp computed from input ikpt and isppol.
!!
!! NOTES
!!
!! PARENTS
!!      newkpt
!!
!! CHILDREN
!!      wffreadskiprec,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine randac(debug,headform1,ikptsp_prev,ikpt,isppol,nband,nkpt,nsppol,wffinp)

 use defs_basis
 use m_wffile
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'randac'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: debug,headform1,ikpt,isppol,nkpt,nsppol
 integer,intent(inout) :: ikptsp_prev
 type(wffile_type),intent(inout) :: wffinp
!arrays
 integer,intent(in) :: nband(nkpt*nsppol)

!Local variables-------------------------------
!scalars
 integer :: ierr,ii,ikptsp,nline,nrec
 character(len=500) :: message

! *************************************************************************

 if (debug>0) then
   write(message, '(a,2i8)' )' randac : ikptsp_prev, ikpt=',ikptsp_prev,ikpt
   call wrtout(std_out,message,'PERS')
 end if

 if(headform1>=40 .or. headform1==0)then
   nline=3 ! npw,npspso,nband record, kg record and eigenvalue record
 else
   nline=2 ! npw,npspso,nband record and eigenvalue record
 end if

!Arrange to skip forward in the data file or else
!backspace as needed to reach the desired records
 ikptsp=ikpt+nkpt*(isppol-1)
 if (ikptsp_prev+1>ikptsp) then
!  Need to backspace nrec records
   nrec=0
   do ii=ikptsp,ikptsp_prev
     nrec=nrec-nband(ii)-nline
   end do
   if (debug>0) then
     write(message, '(a,i8)' )' randac skip back nrec=',-nrec
     call wrtout(std_out,message,'PERS')
   end if
 else if (ikptsp_prev+1<ikptsp) then
!  Need to skip forward nrec records
   nrec=0
   do ii=ikptsp_prev+1,ikptsp-1
!    Need additional 3 again for npw,npspso,nband, kg and eigenvalue records
     nrec=nrec+nband(ii)+nline
   end do
   if (debug>0) then
     write(message, '(a,i8)' )' randac skip forwards nrec=',nrec
     call wrtout(std_out,message,'PERS')
   end if
 else
!  Already pointed at desired record; no skipping
   nrec=0
 end if

!DEBUG
!write(std_out,*)' randac : nrec,wfinp=',nrec,wffinp%unwff
!ENDDEBUG

!Do the skipping (backward, forward, or not)
 call WffReadSkipRec(ierr,nrec,wffinp)

 ikptsp_prev=ikptsp
 if (debug>0) then
   write(message, '(a,i5)' )' randac: updated ikptsp_prev=',ikptsp_prev
   call wrtout(std_out,message,'PERS')
 end if

end subroutine randac
!!***

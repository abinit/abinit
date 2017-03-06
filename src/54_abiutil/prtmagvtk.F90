!{\src2tex{textfont=tt}}
!!****f* ABINIT/prtmagvtk
!! NAME
!!  prtmagvtk
!!
!! FUNCTION
!!  Auxiliary routine for printing out magnetization density in VTK format.
!!
!! COPYRIGHT
!!  Copyright (C) 2017 ABINIT group (SPr)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  argin(sizein)=description
!!
!! OUTPUT
!!  argout(sizeout)=description
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!!  At the moment this routine is needed for development and debugging 
!!  of gs and dfpt calculations with non-collinear spins.
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


subroutine prtmagvtk(nfft,ngfft,nspden,rhor,rprimd)
    
 use defs_basis
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prtmagvtk'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in)  :: nfft,nspden
!arrays
 integer,intent(in)  :: ngfft(18)
 real(dp),intent(in) :: rhor(nfft,nspden),rprimd(3,3)
!real(dp),intent(in) :: xred(3,natom)

!Local variables-------------------------------
!integer ::                                      
!real(dp) ::                                     
!character(len=500) :: msg                   
 
! *************************************************************************

 DBG_ENTER("COLL")
 
! if (option/=1 .and. option/=2 ) then
!   write(msg,'(3a,i0)')&
!&   'The argument option should be 1 or 2,',ch10,&
!&   'however, option=',option
!   MSG_BUG(msg)
! end if
!
! if (sizein<1) then
!   write(msg,'(3a,i0)')&
!&   'The argument sizein should be a positive number,',ch10,&
!&   'however, sizein=',sizein
!   MSG_ERROR(msg)
! end if

 DBG_EXIT("COLL")

 write(238,*) 'entering prtmagvtk, input variables:'
 write(238,*) ' nfft  : ',nspden
 write(238,*) ' nspden: ',nfft
 write(238,*) ' ngfft : ',ngfft(:)
 write(238,*) ' rprimd: '
 write(238,*) '         ',rprimd(1,:)
 write(238,*) '         ',rprimd(2,:)
 write(238,*) '         ',rprimd(3,:)


end subroutine prtmagvtk
!!***

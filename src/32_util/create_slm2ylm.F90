!{\src2tex{textfont=tt}}
!!****f* ABINIT/create_slm2ylm
!! NAME
!! create_slm2ylm
!!
!! FUNCTION
!! For a given angular momentum lcor, compute slm2ylm
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (BA)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  lcor= angular momentum, size of the matrix is 2(2*lcor+1)
!!
!! OUTPUT
!!  slm2ylm(2lcor+1,2lcor+1) = rotation matrix.
!!
!! NOTES
!!  usefull only in ndij==4
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

subroutine create_slm2ylm(lcor,slmtwoylm)

 use defs_basis
 use defs_datatypes
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'create_slm2ylm'
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: lcor
!arrays
 complex(dpc),intent(out) :: slmtwoylm(2*lcor+1,2*lcor+1)

!Local variables ---------------------------------------
!scalars
 integer :: jm,ll,mm,im
 real(dp),parameter :: invsqrt2=one/sqrt2
 real(dp) :: onem
!arrays

! *********************************************************************

 ll=lcor
 slmtwoylm=czero
 do im=1,2*ll+1
   mm=im-ll-1;jm=-mm+ll+1
   onem=dble((-1)**mm)
   if (mm> 0) then
     slmtwoylm(im,im)= cmplx(onem*invsqrt2,zero,kind=dp)
     slmtwoylm(jm,im)= cmplx(invsqrt2,     zero,kind=dp)
   end if
   if (mm==0) then
     slmtwoylm(im,im)=cone
   end if
   if (mm< 0) then
     slmtwoylm(im,im)= cmplx(zero,     invsqrt2,kind=dp)
     slmtwoylm(jm,im)=-cmplx(zero,onem*invsqrt2,kind=dp)
   end if
 end do

end subroutine create_slm2ylm
!!***

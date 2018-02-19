!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkgrid_fft
!! NAME
!!  mkgrid_fft
!!
!! FUNCTION
!!  It sets the grid of fft (or real space) points to be treated.
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2018 ABINIT group (T.Rangel, DC)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      mkcore_paw,mklocl_realspace
!!
!! CHILDREN
!!      xred2xcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine mkgrid_fft(ffti3_local,fftn3_distrib,gridcart,nfft,ngfft,rprimd)
    
 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkgrid_fft'
 use interfaces_41_geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer, intent(in) :: nfft
 integer,intent(in) :: ngfft(18)
 integer, dimension(*), intent(in) :: ffti3_local,fftn3_distrib
 real(dp), dimension(3,nfft), intent(out) :: gridcart
 real(dp),intent(in) :: rprimd(3,3)

!Local variables-------------------------------
 integer :: ind,i1,i2,i3,i3loc,me,nproc
 integer :: n1,n2,n3
 real(dp), dimension(3) :: coord
 real(dp), dimension(3,nfft) :: gridred
!character(len=500) :: message                   ! to be uncommented, if needed
 
! *************************************************************************

 n1    = ngfft(1)
 n2    = ngfft(2)
 n3    = ngfft(3)
 nproc = ngfft(10)
 me    = ngfft(11)

 do i3 = 1, n3, 1
   if(fftn3_distrib(i3) == me) then !MPI
     i3loc=ffti3_local(i3)
     coord(3) = real(i3 - 1, dp) / real(n3, dp)
     do i2 = 1, n2, 1
       coord(2) = real(i2 - 1, dp) / real(n2, dp)
       do i1 = 1, n1, 1
         ind=i1+(i2-1)*n1+(i3loc-1)*n1*n2
         coord(1) = real(i1 - 1, dp) / real(n1, dp)
         gridred(:, ind) = coord(:)
       end do
     end do
   end if
 end do
 call xred2xcart(nfft, rprimd, gridcart, gridred)

 

end subroutine mkgrid_fft
!!***

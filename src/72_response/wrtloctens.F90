!{\src2tex{textfont=tt}}
!!****f* ABINIT/wrtloctens
!! NAME
!! wrtloctens
!!
!! FUNCTION
!! Output of the localisation tensor
!!
!! COPYRIGHT
!! Copyright (C) 1999-2018 ABINIT group (mkv)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! blkflg = flags for each element of the 2DTE (=1 if computed)
!! d2bbb  = band by band decomposition of second order derivatives
!! d2nl   = non-local contributions to the 2DTEs
!! mband  = maximum number of bands
!! mpert  = maximum number of ipert
!! natom  = number of atoms in unit cell
!! prtbbb = if = 1, write the band by band decomposition of the localization tensor
!! rprimd = dimensional primitive translations for real space (bohr)
!! usepaw = flag for PAW
!!
!! OUTPUT
!!  (only writing)
!!
!! TODO
!!  The localization tensor cannot be defined in the metallic case. It should not be computed.
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine wrtloctens(blkflg,d2bbb,d2nl,mband,mpert,natom,prtbbb,rprimd,usepaw)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wrtloctens'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mpert,natom,prtbbb,usepaw
!arrays
 integer,intent(in) :: blkflg(3,mpert,3,mpert)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(inout) :: d2bbb(2,3,3,mpert,mband,mband*prtbbb)
 real(dp),intent(inout) :: d2nl(2,3,mpert,3,mpert)

!Local variables ------------------------------
!scalars
 integer :: flag,iband,idir,idir2,jband
 character(len=500) :: message
!arrays
 real(dp) :: loctenscart(2,3,3)
 real(dp),allocatable :: loctenscart_bbb(:,:,:,:,:)

! *********************************************************************

!This feature is disabled in the PAW case
 if (usepaw==1) return

 if(prtbbb==1)then
   ABI_ALLOCATE(loctenscart_bbb,(2,3,3,mband,mband*prtbbb))
   loctenscart_bbb(:,:,:,:,:)=zero
 end if

!complete missing elements
 flag = 0
 do idir2 = 1,3
   do idir = 1,3

     if (blkflg(idir2,natom+1,idir,natom+1)==0) then
       if (blkflg(idir,natom+1,idir2,natom+1)==0) then
         flag = 1
       else
         d2nl(1,idir2,natom+1,idir,natom+1) = d2nl(1,idir,natom+1,idir2,natom+1)
         d2nl(2,idir2,natom+1,idir,natom+1) =-d2nl(2,idir,natom+1,idir2,natom+1)
         if(prtbbb==1)then
           d2bbb(1,idir2,idir,natom+1,:,:) = d2bbb(1,idir,idir2,natom+1,:,:)
           d2bbb(2,idir2,idir,natom+1,:,:) =-d2bbb(2,idir,idir2,natom+1,:,:)
         end if

       end if
     end if

   end do ! idir=1,3
 end do ! idir2=1,3

!Transform the tensor to cartesian coordinates

 loctenscart(1,:,:) = matmul(rprimd,d2nl(1,:,natom+1,:,natom+1))
 loctenscart(2,:,:) = matmul(rprimd,d2nl(2,:,natom+1,:,natom+1))

 loctenscart(1,:,:) = matmul(loctenscart(1,:,:),transpose(rprimd))
 loctenscart(2,:,:) = matmul(loctenscart(2,:,:),transpose(rprimd))

 loctenscart(:,:,:) = loctenscart(:,:,:)/(two_pi**2)

 if (prtbbb == 1) then

   write(message,'(a,a)')ch10, &
&   ' Band by band decomposition of the localisation tensor (bohr^2)'
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

   do iband = 1,mband
     do jband = 1,mband

       loctenscart_bbb(1,:,:,iband,jband) = matmul(rprimd,d2bbb(1,:,:,natom+1,iband,jband))
       loctenscart_bbb(2,:,:,iband,jband) = matmul(rprimd,d2bbb(2,:,:,natom+1,iband,jband))

       loctenscart_bbb(1,:,:,iband,jband) = matmul(loctenscart_bbb(1,:,:,iband,jband),transpose(rprimd))
       loctenscart_bbb(2,:,:,iband,jband) = matmul(loctenscart_bbb(2,:,:,iband,jband),transpose(rprimd))

       loctenscart_bbb(:,:,:,iband,jband) = loctenscart_bbb(:,:,:,iband,jband)/(two_pi**2)


       write(message,'(a,a,i5,a,i5,a)')ch10, &
&       ' Localisation tensor (bohr^2) for band ',iband,',',jband, &
&       ' in cartesian coordinates'
       call wrtout(std_out,message,'COLL')
       call wrtout(ab_out,message,'COLL')

       write(ab_out,*)'     direction              matrix element'
       write(ab_out,*)'  alpha     beta       real part   imaginary part'
       do idir2 = 1,3
         do idir = 1,3
           write(ab_out,'(5x,i1,8x,i1,3x,2f16.10)')idir2,idir,&
&           loctenscart_bbb(1,idir2,idir,iband,jband),&
&           loctenscart_bbb(2,idir2,idir,iband,jband)
         end do
       end do

     end do !jband
   end do !iband

 end if  !prtbbb

 if (usepaw==0) then
   write(message,'(a,a,a,a)')ch10, &
&   ' Total localisation tensor (bohr^2) in cartesian coordinates',ch10,&
&   '  WARNING : still subject to testing - especially symmetries.'
 else
   write(message,'(a,a,a,a,a,a)')ch10, &
&   ' Total localisation tensor (bohr^2) in cartesian coordinates',ch10,&
&   '  WARNING : probably wrong for PAW (printing for testing purpose)',ch10,&
&   '  WARNING : still subject to testing - especially symmetries.'
 end if
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')

 write(ab_out,*)'     direction              matrix element'
 write(ab_out,*)'  alpha     beta       real part   imaginary part'
 do idir2 = 1,3
   do idir = 1,3
     write(ab_out,'(5x,i1,8x,i1,3x,2f16.10)')idir2,idir,&
&     loctenscart(1,idir2,idir),&
&     loctenscart(2,idir2,idir)
   end do
 end do

 if (flag == 1) then
   write(message,'(6a)')ch10,&
&   ' WARNING : Localization tensor calculation (this does not apply to other properties).',ch10,&
&   '  Not all d/dk perturbations were computed. So the localization tensor in reciprocal space is incomplete,',ch10,&
&   '  and transformation to cartesian coordinates may be wrong.'
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')
 end if

 if (prtbbb == 1) then
   ABI_DEALLOCATE(loctenscart_bbb)
 end if

end subroutine wrtloctens
!!***

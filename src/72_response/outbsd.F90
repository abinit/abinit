!{\src2tex{textfont=tt}}
!!****f* ABINIT/outbsd
!! NAME
!! outbsd
!!
!! FUNCTION
!! output bsd file for one perturbation (used for elphon calculations in anaddb)
!!
!! COPYRIGHT
!! Copyright (C) 1999-2017 ABINIT group (XG, DRH, MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.
!!
!! INPUTS
!!  bdeigrf=number of bands for which the derivatives of the eigenvalues have been computed
!!  dtset = dataset variable for run flags
!!  eig2nkq= second ordre eigenvalue (or electron lifetime) that must be printed out 
!!  mpert= maximum number of perturbations
!!  nkpt_rbz= number of k-points for perturbation
!!  unitout= writting unit of file
!!
!! OUTPUTS
!!  to file
!!
!! PARENTS
!!      dfpt_looppert,eig2tot
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine outbsd(bdeigrf,dtset,eig2nkq,mpert,nkpt_rbz,unitout)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'outbsd'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: bdeigrf,mpert,nkpt_rbz,unitout
 type(dataset_type),intent(in) :: dtset
!arrays
 real(dp),intent(in) :: eig2nkq(2,dtset%mband*dtset%nsppol,nkpt_rbz,3,mpert,3,mpert)

!Local variables -------------------------
!scalars
 integer :: bandtot_index,iband,idir1,idir2,ikpt,ipert1,ipert2,isppol

! *********************************************************************

!DEBUG
!write(std_out,*)' outbsd : enter'
!write(std_out,*)' eig2nkq(1,1,1,1,1,1)=',eig2nkq(1,1,1,1,1,1,1)
!ENDDEBUG

!output information in this file
 write(unitout,*)
 write(unitout,'(a,i8)') ' 2nd eigenvalue derivatives   - # elements :', 9*dtset%natom**2
 write(unitout,'(a,3es16.8,a)') ' qpt', dtset%qptn(:), ' 1.0'

!output RF eigenvalues

 do ikpt=1,nkpt_rbz
!  bandtot_index differs from zero only in the spin-polarized case
   bandtot_index=0
   write (unitout,'(a,3es16.8)') ' K-point:', dtset%kptns(:,ikpt)
   do isppol=1,dtset%nsppol
     do iband=1,bdeigrf
       write (unitout,'(a,i5)') ' Band:', iband+bandtot_index
!      write (unitout,*) 'ipert1     ','idir1     ','ipert2     ','idir2    ','Real    ','Im    '
       do ipert2=1,mpert
         do idir2=1,3
           do ipert1=1,mpert
             do idir1=1,3
               write (unitout,'(4i4,2d22.14)') idir1,ipert1,idir2,ipert2,&
&               eig2nkq(1,iband+bandtot_index,ikpt,idir1,ipert1,idir2,ipert2),&
&               eig2nkq(2,iband+bandtot_index,ikpt,idir1,ipert1,idir2,ipert2)
             end do !idir2
           end do !ipert2
         end do !idir1
       end do !ipert1
     end do !iband
     bandtot_index = bandtot_index + dtset%mband
   end do !isppol
 end do !ikpt

!close bsd file
 close (unitout)

end subroutine outbsd
!!***

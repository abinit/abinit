!{\src2tex{textfont=tt}}
!!****f* ABINIT/d2c_wtq
!! NAME
!!  d2c_wtq
!!
!! FUNCTION
!!  This routine calculates the integration weights on a coarse k-grid 
!!   using the integration weights from a denser k-gird. The weights of 
!!   the extra k points that being shared are evenly distributed.
!!
!! COPYRIGHT
!!  Copyright (C) 2012-2018 ABINIT group (BXU)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  elph_ds%k_fine%nkpt = number of fine q-points
!!  elph_ds%k_fine%wtq = integration weights of the fine q-grid
!!  elph_ds%k_phon%nkpt = number of coarse q-points
!!
!! OUTPUT
!!  elph_ds%k_phon%wtq = integration weights of the coarse k-grid
!!
!! PARENTS
!!      mka2f,mka2f_tr
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine d2c_wtq(elph_ds)
    
 use defs_basis
 use defs_elphon
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'd2c_wtq'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(elph_type),intent(inout) :: elph_ds

!Local variables-------------------------------
 integer :: ii, jj, kk 
 integer :: ikpt, jkpt, kkpt
 integer :: iikpt, jjkpt, kkkpt
 integer :: ikpt_fine, ikpt_phon
 integer :: nkpt_fine1, nkpt_phon1
 integer :: nkpt_fine2, nkpt_phon2
 integer :: nkpt_fine3, nkpt_phon3
 integer :: nscale1, nscale2, nscale3
 
! *************************************************************************
 nkpt_phon1 = elph_ds%kptrlatt(1,1)
 nkpt_phon2 = elph_ds%kptrlatt(2,2)
 nkpt_phon3 = elph_ds%kptrlatt(3,3)
 nkpt_fine1 = elph_ds%kptrlatt_fine(1,1)
 nkpt_fine2 = elph_ds%kptrlatt_fine(2,2)
 nkpt_fine3 = elph_ds%kptrlatt_fine(3,3)
 nscale1 = dble(nkpt_fine1/nkpt_phon1)
 nscale2 = dble(nkpt_fine2/nkpt_phon2)
 nscale3 = dble(nkpt_fine3/nkpt_phon3)
 if (abs(INT(nscale1)-nscale1) > 0.01) then
   MSG_ERROR('The denser k-gird MUST be multiples of the phon k-grid')
 end if
 if (abs(INT(nscale2)-nscale2) > 0.01) then
   MSG_ERROR('The denser k-gird MUST be multiples of the phon k-grid')
 end if
 if (abs(INT(nscale3)-nscale3) > 0.01) then
   MSG_ERROR('The denser k-gird MUST be multiples of the phon k-grid')
 end if
 nscale1 = INT(nscale1)
 nscale2 = INT(nscale2)
 nscale3 = INT(nscale3)
 
!bxu, get wtq of coarse grid from fine grid
 elph_ds%k_phon%wtq = zero

 do ikpt = 1, nkpt_phon1
   do jkpt = 1, nkpt_phon2
     do kkpt = 1, nkpt_phon3
       ikpt_phon = kkpt + (jkpt-1)*nkpt_phon3 + (ikpt-1)*nkpt_phon2*nkpt_phon3
!      inside the paralellepipe
       do ii = -((nscale1+1)/2-1), ((nscale1+1)/2-1)
         do jj = -((nscale2+1)/2-1), ((nscale2+1)/2-1)
           do kk = -((nscale3+1)/2-1), ((nscale3+1)/2-1)
             iikpt = 1 + (ikpt-1)*nscale1 + ii
             jjkpt = 1 + (jkpt-1)*nscale2 + jj
             kkkpt = 1 + (kkpt-1)*nscale3 + kk
             if (iikpt .le. 0) iikpt = iikpt + nkpt_fine1
             if (jjkpt .le. 0) jjkpt = jjkpt + nkpt_fine2
             if (kkkpt .le. 0) kkkpt = kkkpt + nkpt_fine3
             if (iikpt .gt. nkpt_fine1) iikpt = iikpt - nkpt_fine1
             if (jjkpt .gt. nkpt_fine2) jjkpt = jjkpt - nkpt_fine2
             if (kkkpt .gt. nkpt_fine3) kkkpt = kkkpt - nkpt_fine3
             ikpt_fine = kkkpt + (jjkpt-1)*nkpt_fine3 + (iikpt-1)*nkpt_fine2*nkpt_fine3
             elph_ds%k_phon%wtq(:,ikpt_phon,:)=elph_ds%k_phon%wtq(:,ikpt_phon,:)+&
&             elph_ds%k_fine%wtq(:,ikpt_fine,:)
           end do
         end do
       end do
!      on the 6 faces
       if (MOD(nscale3,2) == 0) then ! when nscale3 is an even number
         do ii = -((nscale1+1)/2-1), ((nscale1+1)/2-1)
           do jj = -((nscale2+1)/2-1), ((nscale2+1)/2-1)
             iikpt = 1 + (ikpt-1)*nscale1 + ii
             jjkpt = 1 + (jkpt-1)*nscale2 + jj
             if (iikpt .le. 0) iikpt = iikpt + nkpt_fine1
             if (jjkpt .le. 0) jjkpt = jjkpt + nkpt_fine2
             if (iikpt .gt. nkpt_fine1) iikpt = iikpt - nkpt_fine1
             if (jjkpt .gt. nkpt_fine2) jjkpt = jjkpt - nkpt_fine2

             kkkpt = 1 + (kkpt-1)*nscale3 + nscale3/2
             if (kkkpt .le. 0) kkkpt = kkkpt + nkpt_fine3
             if (kkkpt .gt. nkpt_fine3) kkkpt = kkkpt - nkpt_fine3
             ikpt_fine = kkkpt + (jjkpt-1)*nkpt_fine3 + (iikpt-1)*nkpt_fine2*nkpt_fine3
             elph_ds%k_phon%wtq(:,ikpt_phon,:)=elph_ds%k_phon%wtq(:,ikpt_phon,:)+&
&             0.5_dp*elph_ds%k_fine%wtq(:,ikpt_fine,:)

             kkkpt = 1 + (kkpt-1)*nscale3 - nscale3/2
             if (kkkpt .le. 0) kkkpt = kkkpt + nkpt_fine3
             if (kkkpt .gt. nkpt_fine3) kkkpt = kkkpt - nkpt_fine3
             ikpt_fine = kkkpt + (jjkpt-1)*nkpt_fine3 + (iikpt-1)*nkpt_fine2*nkpt_fine3
             elph_ds%k_phon%wtq(:,ikpt_phon,:)=elph_ds%k_phon%wtq(:,ikpt_phon,:)+&
&             0.5_dp*elph_ds%k_fine%wtq(:,ikpt_fine,:)
           end do
         end do
       end if 
       if (MOD(nscale2,2) == 0) then ! when nscale2 is an even number
         do ii = -((nscale1+1)/2-1), ((nscale1+1)/2-1)
           do kk = -((nscale3+1)/2-1), ((nscale3+1)/2-1)
             iikpt = 1 + (ikpt-1)*nscale1 + ii
             kkkpt = 1 + (kkpt-1)*nscale3 + kk
             if (iikpt .le. 0) iikpt = iikpt + nkpt_fine1
             if (kkkpt .le. 0) kkkpt = kkkpt + nkpt_fine3
             if (iikpt .gt. nkpt_fine1) iikpt = iikpt - nkpt_fine1
             if (kkkpt .gt. nkpt_fine3) kkkpt = kkkpt - nkpt_fine3

             jjkpt = 1 + (jkpt-1)*nscale2 + nscale2/2
             if (jjkpt .le. 0) jjkpt = jjkpt + nkpt_fine2
             if (jjkpt .gt. nkpt_fine2) jjkpt = jjkpt - nkpt_fine2
             ikpt_fine = kkkpt + (jjkpt-1)*nkpt_fine3 + (iikpt-1)*nkpt_fine2*nkpt_fine3
             elph_ds%k_phon%wtq(:,ikpt_phon,:)=elph_ds%k_phon%wtq(:,ikpt_phon,:)+&
&             0.5_dp*elph_ds%k_fine%wtq(:,ikpt_fine,:)

             jjkpt = 1 + (jkpt-1)*nscale2 - nscale2/2
             if (jjkpt .le. 0) jjkpt = jjkpt + nkpt_fine2
             if (jjkpt .gt. nkpt_fine2) jjkpt = jjkpt - nkpt_fine2
             ikpt_fine = kkkpt + (jjkpt-1)*nkpt_fine3 + (iikpt-1)*nkpt_fine2*nkpt_fine3
             elph_ds%k_phon%wtq(:,ikpt_phon,:)=elph_ds%k_phon%wtq(:,ikpt_phon,:)+&
&             0.5_dp*elph_ds%k_fine%wtq(:,ikpt_fine,:)
           end do
         end do
       end if 
       if (MOD(nscale1,2) == 0) then ! when nscale1 is an even number
         do kk = -((nscale3+1)/2-1), ((nscale3+1)/2-1)
           do jj = -((nscale2+1)/2-1), ((nscale2+1)/2-1)
             kkkpt = 1 + (kkpt-1)*nscale3 + kk
             jjkpt = 1 + (jkpt-1)*nscale2 + jj
             if (kkkpt .le. 0) kkkpt = kkkpt + nkpt_fine3
             if (jjkpt .le. 0) jjkpt = jjkpt + nkpt_fine2
             if (kkkpt .gt. nkpt_fine3) kkkpt = kkkpt - nkpt_fine3
             if (jjkpt .gt. nkpt_fine2) jjkpt = jjkpt - nkpt_fine2

             iikpt = 1 + (ikpt-1)*nscale1 + nscale1/2
             if (iikpt .le. 0) iikpt = iikpt + nkpt_fine1
             if (iikpt .gt. nkpt_fine1) iikpt = iikpt - nkpt_fine1
             ikpt_fine = kkkpt + (jjkpt-1)*nkpt_fine3 + (iikpt-1)*nkpt_fine2*nkpt_fine3
             elph_ds%k_phon%wtq(:,ikpt_phon,:)=elph_ds%k_phon%wtq(:,ikpt_phon,:)+&
&             0.5_dp*elph_ds%k_fine%wtq(:,ikpt_fine,:)

             iikpt = 1 + (ikpt-1)*nscale1 - nscale1/2
             if (iikpt .le. 0) iikpt = iikpt + nkpt_fine1
             if (iikpt .gt. nkpt_fine1) iikpt = iikpt - nkpt_fine1
             ikpt_fine = kkkpt + (jjkpt-1)*nkpt_fine3 + (iikpt-1)*nkpt_fine2*nkpt_fine3
             elph_ds%k_phon%wtq(:,ikpt_phon,:)=elph_ds%k_phon%wtq(:,ikpt_phon,:)+&
&             0.5_dp*elph_ds%k_fine%wtq(:,ikpt_fine,:)
           end do
         end do
!        on the 12 sides
       end if 
       if (MOD(nscale2,2) == 0 .and. MOD(nscale3,2) == 0) then
         do ii = -((nscale1+1)/2-1), ((nscale1+1)/2-1)
           iikpt = 1 + (ikpt-1)*nscale1 + ii
           if (iikpt .le. 0) iikpt = iikpt + nkpt_fine1
           if (iikpt .gt. nkpt_fine1) iikpt = iikpt - nkpt_fine1

           jjkpt = 1 + (jkpt-1)*nscale2 + nscale2/2
           kkkpt = 1 + (kkpt-1)*nscale3 + nscale3/2
           if (jjkpt .le. 0) jjkpt = jjkpt + nkpt_fine2
           if (jjkpt .gt. nkpt_fine2) jjkpt = jjkpt - nkpt_fine2
           if (kkkpt .le. 0) kkkpt = kkkpt + nkpt_fine3
           if (kkkpt .gt. nkpt_fine3) kkkpt = kkkpt - nkpt_fine3
           ikpt_fine = kkkpt + (jjkpt-1)*nkpt_fine3 + (iikpt-1)*nkpt_fine2*nkpt_fine3
           elph_ds%k_phon%wtq(:,ikpt_phon,:)=elph_ds%k_phon%wtq(:,ikpt_phon,:)+&
&           0.25_dp*elph_ds%k_fine%wtq(:,ikpt_fine,:)

           jjkpt = 1 + (jkpt-1)*nscale2 + nscale2/2
           kkkpt = 1 + (kkpt-1)*nscale3 - nscale3/2
           if (jjkpt .le. 0) jjkpt = jjkpt + nkpt_fine2
           if (jjkpt .gt. nkpt_fine2) jjkpt = jjkpt - nkpt_fine2
           if (kkkpt .le. 0) kkkpt = kkkpt + nkpt_fine3
           if (kkkpt .gt. nkpt_fine3) kkkpt = kkkpt - nkpt_fine3
           ikpt_fine = kkkpt + (jjkpt-1)*nkpt_fine3 + (iikpt-1)*nkpt_fine2*nkpt_fine3
           elph_ds%k_phon%wtq(:,ikpt_phon,:)=elph_ds%k_phon%wtq(:,ikpt_phon,:)+&
&           0.25_dp*elph_ds%k_fine%wtq(:,ikpt_fine,:)

           jjkpt = 1 + (jkpt-1)*nscale2 - nscale2/2
           kkkpt = 1 + (kkpt-1)*nscale3 + nscale3/2
           if (jjkpt .le. 0) jjkpt = jjkpt + nkpt_fine2
           if (jjkpt .gt. nkpt_fine2) jjkpt = jjkpt - nkpt_fine2
           if (kkkpt .le. 0) kkkpt = kkkpt + nkpt_fine3
           if (kkkpt .gt. nkpt_fine3) kkkpt = kkkpt - nkpt_fine3
           ikpt_fine = kkkpt + (jjkpt-1)*nkpt_fine3 + (iikpt-1)*nkpt_fine2*nkpt_fine3
           elph_ds%k_phon%wtq(:,ikpt_phon,:)=elph_ds%k_phon%wtq(:,ikpt_phon,:)+&
&           0.25_dp*elph_ds%k_fine%wtq(:,ikpt_fine,:)

           jjkpt = 1 + (jkpt-1)*nscale2 - nscale2/2
           kkkpt = 1 + (kkpt-1)*nscale3 - nscale3/2
           if (jjkpt .le. 0) jjkpt = jjkpt + nkpt_fine2
           if (jjkpt .gt. nkpt_fine2) jjkpt = jjkpt - nkpt_fine2
           if (kkkpt .le. 0) kkkpt = kkkpt + nkpt_fine3
           if (kkkpt .gt. nkpt_fine3) kkkpt = kkkpt - nkpt_fine3
           ikpt_fine = kkkpt + (jjkpt-1)*nkpt_fine3 + (iikpt-1)*nkpt_fine2*nkpt_fine3
           elph_ds%k_phon%wtq(:,ikpt_phon,:)=elph_ds%k_phon%wtq(:,ikpt_phon,:)+&
&           0.25_dp*elph_ds%k_fine%wtq(:,ikpt_fine,:)
         end do
       end if 
       if (MOD(nscale1,2) == 0 .and. MOD(nscale3,2) == 0) then
         do jj = -((nscale2+1)/2-1), ((nscale2+1)/2-1)
           jjkpt = 1 + (jkpt-1)*nscale2 + jj
           if (jjkpt .le. 0) jjkpt = jjkpt + nkpt_fine2
           if (jjkpt .gt. nkpt_fine2) jjkpt = jjkpt - nkpt_fine2

           iikpt = 1 + (ikpt-1)*nscale1 + nscale1/2
           kkkpt = 1 + (kkpt-1)*nscale3 + nscale3/2
           if (iikpt .le. 0) iikpt = iikpt + nkpt_fine1
           if (iikpt .gt. nkpt_fine1) iikpt = iikpt - nkpt_fine1
           if (kkkpt .le. 0) kkkpt = kkkpt + nkpt_fine3
           if (kkkpt .gt. nkpt_fine3) kkkpt = kkkpt - nkpt_fine3
           ikpt_fine = kkkpt + (jjkpt-1)*nkpt_fine3 + (iikpt-1)*nkpt_fine2*nkpt_fine3
           elph_ds%k_phon%wtq(:,ikpt_phon,:)=elph_ds%k_phon%wtq(:,ikpt_phon,:)+&
&           0.25_dp*elph_ds%k_fine%wtq(:,ikpt_fine,:)

           iikpt = 1 + (ikpt-1)*nscale1 + nscale1/2
           kkkpt = 1 + (kkpt-1)*nscale3 - nscale3/2
           if (iikpt .le. 0) iikpt = iikpt + nkpt_fine1
           if (iikpt .gt. nkpt_fine1) iikpt = iikpt - nkpt_fine1
           if (kkkpt .le. 0) kkkpt = kkkpt + nkpt_fine3
           if (kkkpt .gt. nkpt_fine3) kkkpt = kkkpt - nkpt_fine3
           ikpt_fine = kkkpt + (jjkpt-1)*nkpt_fine3 + (iikpt-1)*nkpt_fine2*nkpt_fine3
           elph_ds%k_phon%wtq(:,ikpt_phon,:)=elph_ds%k_phon%wtq(:,ikpt_phon,:)+&
&           0.25_dp*elph_ds%k_fine%wtq(:,ikpt_fine,:)

           iikpt = 1 + (ikpt-1)*nscale1 - nscale1/2
           kkkpt = 1 + (kkpt-1)*nscale3 + nscale3/2
           if (iikpt .le. 0) iikpt = iikpt + nkpt_fine1
           if (iikpt .gt. nkpt_fine1) iikpt = iikpt - nkpt_fine1
           if (kkkpt .le. 0) kkkpt = kkkpt + nkpt_fine3
           if (kkkpt .gt. nkpt_fine3) kkkpt = kkkpt - nkpt_fine3
           ikpt_fine = kkkpt + (jjkpt-1)*nkpt_fine3 + (iikpt-1)*nkpt_fine2*nkpt_fine3
           elph_ds%k_phon%wtq(:,ikpt_phon,:)=elph_ds%k_phon%wtq(:,ikpt_phon,:)+&
&           0.25_dp*elph_ds%k_fine%wtq(:,ikpt_fine,:)

           iikpt = 1 + (ikpt-1)*nscale1 - nscale1/2
           kkkpt = 1 + (kkpt-1)*nscale3 - nscale3/2
           if (iikpt .le. 0) iikpt = iikpt + nkpt_fine1
           if (iikpt .gt. nkpt_fine1) iikpt = iikpt - nkpt_fine1
           if (kkkpt .le. 0) kkkpt = kkkpt + nkpt_fine3
           if (kkkpt .gt. nkpt_fine3) kkkpt = kkkpt - nkpt_fine3
           ikpt_fine = kkkpt + (jjkpt-1)*nkpt_fine3 + (iikpt-1)*nkpt_fine2*nkpt_fine3
           elph_ds%k_phon%wtq(:,ikpt_phon,:)=elph_ds%k_phon%wtq(:,ikpt_phon,:)+&
&           0.25_dp*elph_ds%k_fine%wtq(:,ikpt_fine,:)
         end do
       end if 
       if (MOD(nscale2,2) == 0 .and. MOD(nscale1,2) == 0) then
         do kk = -((nscale3+1)/2-1), ((nscale3+1)/2-1)
           kkkpt = 1 + (kkpt-1)*nscale3 + kk
           if (kkkpt .le. 0) kkkpt = kkkpt + nkpt_fine3
           if (kkkpt .gt. nkpt_fine3) kkkpt = kkkpt - nkpt_fine3

           jjkpt = 1 + (jkpt-1)*nscale2 + nscale2/2
           iikpt = 1 + (ikpt-1)*nscale1 + nscale1/2
           if (jjkpt .le. 0) jjkpt = jjkpt + nkpt_fine2
           if (jjkpt .gt. nkpt_fine2) jjkpt = jjkpt - nkpt_fine2
           if (iikpt .le. 0) iikpt = iikpt + nkpt_fine1
           if (iikpt .gt. nkpt_fine1) iikpt = iikpt - nkpt_fine1
           ikpt_fine = kkkpt + (jjkpt-1)*nkpt_fine3 + (iikpt-1)*nkpt_fine2*nkpt_fine3
           elph_ds%k_phon%wtq(:,ikpt_phon,:)=elph_ds%k_phon%wtq(:,ikpt_phon,:)+&
&           0.25_dp*elph_ds%k_fine%wtq(:,ikpt_fine,:)

           jjkpt = 1 + (jkpt-1)*nscale2 + nscale2/2
           iikpt = 1 + (ikpt-1)*nscale1 - nscale1/2
           if (jjkpt .le. 0) jjkpt = jjkpt + nkpt_fine2
           if (jjkpt .gt. nkpt_fine2) jjkpt = jjkpt - nkpt_fine2
           if (iikpt .le. 0) iikpt = iikpt + nkpt_fine1
           if (iikpt .gt. nkpt_fine1) iikpt = iikpt - nkpt_fine1
           ikpt_fine = kkkpt + (jjkpt-1)*nkpt_fine3 + (iikpt-1)*nkpt_fine2*nkpt_fine3
           elph_ds%k_phon%wtq(:,ikpt_phon,:)=elph_ds%k_phon%wtq(:,ikpt_phon,:)+&
&           0.25_dp*elph_ds%k_fine%wtq(:,ikpt_fine,:)

           jjkpt = 1 + (jkpt-1)*nscale2 - nscale2/2
           iikpt = 1 + (ikpt-1)*nscale1 + nscale1/2
           if (jjkpt .le. 0) jjkpt = jjkpt + nkpt_fine2
           if (jjkpt .gt. nkpt_fine2) jjkpt = jjkpt - nkpt_fine2
           if (iikpt .le. 0) iikpt = iikpt + nkpt_fine1
           if (iikpt .gt. nkpt_fine1) iikpt = iikpt - nkpt_fine1
           ikpt_fine = kkkpt + (jjkpt-1)*nkpt_fine3 + (iikpt-1)*nkpt_fine2*nkpt_fine3
           elph_ds%k_phon%wtq(:,ikpt_phon,:)=elph_ds%k_phon%wtq(:,ikpt_phon,:)+&
&           0.25_dp*elph_ds%k_fine%wtq(:,ikpt_fine,:)

           jjkpt = 1 + (jkpt-1)*nscale2 - nscale2/2
           iikpt = 1 + (ikpt-1)*nscale1 - nscale1/2
           if (jjkpt .le. 0) jjkpt = jjkpt + nkpt_fine2
           if (jjkpt .gt. nkpt_fine2) jjkpt = jjkpt - nkpt_fine2
           if (iikpt .le. 0) iikpt = iikpt + nkpt_fine1
           if (iikpt .gt. nkpt_fine1) iikpt = iikpt - nkpt_fine1
           ikpt_fine = kkkpt + (jjkpt-1)*nkpt_fine3 + (iikpt-1)*nkpt_fine2*nkpt_fine3
           elph_ds%k_phon%wtq(:,ikpt_phon,:)=elph_ds%k_phon%wtq(:,ikpt_phon,:)+&
&           0.25_dp*elph_ds%k_fine%wtq(:,ikpt_fine,:)
         end do
!        on the 8 corners
       end if 
       if (MOD(nscale1,2) == 0 .and. MOD(nscale2,2) == 0 .and. MOD(nscale3,2) == 0) then
         iikpt = 1 + (ikpt-1)*nscale1 + nscale1/2
         jjkpt = 1 + (jkpt-1)*nscale2 + nscale2/2
         kkkpt = 1 + (kkpt-1)*nscale3 + nscale3/2
         if (iikpt .le. 0) iikpt = iikpt + nkpt_fine1
         if (jjkpt .le. 0) jjkpt = jjkpt + nkpt_fine2
         if (kkkpt .le. 0) kkkpt = kkkpt + nkpt_fine3
         if (iikpt .gt. nkpt_fine1) iikpt = iikpt - nkpt_fine1
         if (jjkpt .gt. nkpt_fine2) jjkpt = jjkpt - nkpt_fine2
         if (kkkpt .gt. nkpt_fine3) kkkpt = kkkpt - nkpt_fine3
         ikpt_fine = kkkpt + (jjkpt-1)*nkpt_fine3 + (iikpt-1)*nkpt_fine2*nkpt_fine3
         elph_ds%k_phon%wtq(:,ikpt_phon,:)=elph_ds%k_phon%wtq(:,ikpt_phon,:)+&
&         0.125_dp*elph_ds%k_fine%wtq(:,ikpt_fine,:)

         iikpt = 1 + (ikpt-1)*nscale1 + nscale1/2
         jjkpt = 1 + (jkpt-1)*nscale2 + nscale2/2
         kkkpt = 1 + (kkpt-1)*nscale3 - nscale3/2
         if (iikpt .le. 0) iikpt = iikpt + nkpt_fine1
         if (jjkpt .le. 0) jjkpt = jjkpt + nkpt_fine2
         if (kkkpt .le. 0) kkkpt = kkkpt + nkpt_fine3
         if (iikpt .gt. nkpt_fine1) iikpt = iikpt - nkpt_fine1
         if (jjkpt .gt. nkpt_fine2) jjkpt = jjkpt - nkpt_fine2
         if (kkkpt .gt. nkpt_fine3) kkkpt = kkkpt - nkpt_fine3
         ikpt_fine = kkkpt + (jjkpt-1)*nkpt_fine3 + (iikpt-1)*nkpt_fine2*nkpt_fine3
         elph_ds%k_phon%wtq(:,ikpt_phon,:)=elph_ds%k_phon%wtq(:,ikpt_phon,:)+&
&         0.125_dp*elph_ds%k_fine%wtq(:,ikpt_fine,:)

         iikpt = 1 + (ikpt-1)*nscale1 + nscale1/2
         jjkpt = 1 + (jkpt-1)*nscale2 - nscale2/2
         kkkpt = 1 + (kkpt-1)*nscale3 + nscale3/2
         if (iikpt .le. 0) iikpt = iikpt + nkpt_fine1
         if (jjkpt .le. 0) jjkpt = jjkpt + nkpt_fine2
         if (kkkpt .le. 0) kkkpt = kkkpt + nkpt_fine3
         if (iikpt .gt. nkpt_fine1) iikpt = iikpt - nkpt_fine1
         if (jjkpt .gt. nkpt_fine2) jjkpt = jjkpt - nkpt_fine2
         if (kkkpt .gt. nkpt_fine3) kkkpt = kkkpt - nkpt_fine3
         ikpt_fine = kkkpt + (jjkpt-1)*nkpt_fine3 + (iikpt-1)*nkpt_fine2*nkpt_fine3
         elph_ds%k_phon%wtq(:,ikpt_phon,:)=elph_ds%k_phon%wtq(:,ikpt_phon,:)+&
&         0.125_dp*elph_ds%k_fine%wtq(:,ikpt_fine,:)

         iikpt = 1 + (ikpt-1)*nscale1 + nscale1/2
         jjkpt = 1 + (jkpt-1)*nscale2 - nscale2/2
         kkkpt = 1 + (kkpt-1)*nscale3 - nscale3/2
         if (iikpt .le. 0) iikpt = iikpt + nkpt_fine1
         if (jjkpt .le. 0) jjkpt = jjkpt + nkpt_fine2
         if (kkkpt .le. 0) kkkpt = kkkpt + nkpt_fine3
         if (iikpt .gt. nkpt_fine1) iikpt = iikpt - nkpt_fine1
         if (jjkpt .gt. nkpt_fine2) jjkpt = jjkpt - nkpt_fine2
         if (kkkpt .gt. nkpt_fine3) kkkpt = kkkpt - nkpt_fine3
         ikpt_fine = kkkpt + (jjkpt-1)*nkpt_fine3 + (iikpt-1)*nkpt_fine2*nkpt_fine3
         elph_ds%k_phon%wtq(:,ikpt_phon,:)=elph_ds%k_phon%wtq(:,ikpt_phon,:)+&
&         0.125_dp*elph_ds%k_fine%wtq(:,ikpt_fine,:)

         iikpt = 1 + (ikpt-1)*nscale1 - nscale1/2
         jjkpt = 1 + (jkpt-1)*nscale2 + nscale2/2
         kkkpt = 1 + (kkpt-1)*nscale3 + nscale3/2
         if (iikpt .le. 0) iikpt = iikpt + nkpt_fine1
         if (jjkpt .le. 0) jjkpt = jjkpt + nkpt_fine2
         if (kkkpt .le. 0) kkkpt = kkkpt + nkpt_fine3
         if (iikpt .gt. nkpt_fine1) iikpt = iikpt - nkpt_fine1
         if (jjkpt .gt. nkpt_fine2) jjkpt = jjkpt - nkpt_fine2
         if (kkkpt .gt. nkpt_fine3) kkkpt = kkkpt - nkpt_fine3
         ikpt_fine = kkkpt + (jjkpt-1)*nkpt_fine3 + (iikpt-1)*nkpt_fine2*nkpt_fine3
         elph_ds%k_phon%wtq(:,ikpt_phon,:)=elph_ds%k_phon%wtq(:,ikpt_phon,:)+&
&         0.125_dp*elph_ds%k_fine%wtq(:,ikpt_fine,:)

         iikpt = 1 + (ikpt-1)*nscale1 - nscale1/2
         jjkpt = 1 + (jkpt-1)*nscale2 + nscale2/2
         kkkpt = 1 + (kkpt-1)*nscale3 - nscale3/2
         if (iikpt .le. 0) iikpt = iikpt + nkpt_fine1
         if (jjkpt .le. 0) jjkpt = jjkpt + nkpt_fine2
         if (kkkpt .le. 0) kkkpt = kkkpt + nkpt_fine3
         if (iikpt .gt. nkpt_fine1) iikpt = iikpt - nkpt_fine1
         if (jjkpt .gt. nkpt_fine2) jjkpt = jjkpt - nkpt_fine2
         if (kkkpt .gt. nkpt_fine3) kkkpt = kkkpt - nkpt_fine3
         ikpt_fine = kkkpt + (jjkpt-1)*nkpt_fine3 + (iikpt-1)*nkpt_fine2*nkpt_fine3
         elph_ds%k_phon%wtq(:,ikpt_phon,:)=elph_ds%k_phon%wtq(:,ikpt_phon,:)+&
&         0.125_dp*elph_ds%k_fine%wtq(:,ikpt_fine,:)

         iikpt = 1 + (ikpt-1)*nscale1 - nscale1/2
         jjkpt = 1 + (jkpt-1)*nscale2 - nscale2/2
         kkkpt = 1 + (kkpt-1)*nscale3 + nscale3/2
         if (iikpt .le. 0) iikpt = iikpt + nkpt_fine1
         if (jjkpt .le. 0) jjkpt = jjkpt + nkpt_fine2
         if (kkkpt .le. 0) kkkpt = kkkpt + nkpt_fine3
         if (iikpt .gt. nkpt_fine1) iikpt = iikpt - nkpt_fine1
         if (jjkpt .gt. nkpt_fine2) jjkpt = jjkpt - nkpt_fine2
         if (kkkpt .gt. nkpt_fine3) kkkpt = kkkpt - nkpt_fine3
         ikpt_fine = kkkpt + (jjkpt-1)*nkpt_fine3 + (iikpt-1)*nkpt_fine2*nkpt_fine3
         elph_ds%k_phon%wtq(:,ikpt_phon,:)=elph_ds%k_phon%wtq(:,ikpt_phon,:)+&
&         0.125_dp*elph_ds%k_fine%wtq(:,ikpt_fine,:)

         iikpt = 1 + (ikpt-1)*nscale1 - nscale1/2
         jjkpt = 1 + (jkpt-1)*nscale2 - nscale2/2
         kkkpt = 1 + (kkpt-1)*nscale3 - nscale3/2
         if (iikpt .le. 0) iikpt = iikpt + nkpt_fine1
         if (jjkpt .le. 0) jjkpt = jjkpt + nkpt_fine2
         if (kkkpt .le. 0) kkkpt = kkkpt + nkpt_fine3
         if (iikpt .gt. nkpt_fine1) iikpt = iikpt - nkpt_fine1
         if (jjkpt .gt. nkpt_fine2) jjkpt = jjkpt - nkpt_fine2
         if (kkkpt .gt. nkpt_fine3) kkkpt = kkkpt - nkpt_fine3
         ikpt_fine = kkkpt + (jjkpt-1)*nkpt_fine3 + (iikpt-1)*nkpt_fine2*nkpt_fine3
         elph_ds%k_phon%wtq(:,ikpt_phon,:)=elph_ds%k_phon%wtq(:,ikpt_phon,:)+&
&         0.125_dp*elph_ds%k_fine%wtq(:,ikpt_fine,:)
       end if
     end do
   end do
 end do

!bxu, divide by nscale^3 to be consistent with the normalization of kpt_phon
 elph_ds%k_phon%wtq = elph_ds%k_phon%wtq/nscale1/nscale2/nscale3

end subroutine d2c_wtq
!!***

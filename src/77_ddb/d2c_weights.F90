!{\src2tex{textfont=tt}}
!!****f* ABINIT/d2c_weights
!! NAME
!!  d2c_weights
!!
!! FUNCTION
!!  This routine calculates the integration weights on a coarse k-grid 
!!   using the integration weights from a denser k-gird. The weights of 
!!   the extra k points that being shared are evenly distributed. It also
!!   condenses the velocity*wtk and velcity^2*wtk.
!!
!! COPYRIGHT
!!  Copyright (C) 2012-2018 ABINIT group (BXU)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  elph_ds%k_fine%nkpt = number of fine FS k-points
!!  elph_ds%k_fine%wtk = integration weights of the fine FS k-grid
!!  elph_ds%k_phon%nkpt = number of coarse FS k-points
!!  elph_tr_ds%el_veloc = electronic velocities from the fine k-grid
!!
!! OUTPUT
!!  elph_ds%k_phon%wtk = integration weights of the coarse FS k-grid
!!  elph_ds%k_phon%velocwtk = velocity time integration weights of the coarse FS k-grid
!!  elph_ds%k_phon%vvelocwtk = velocity^2 time integration weights of the coarse FS k-grid
!!
!! PARENTS
!!      elphon,get_nv_fs_en
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine d2c_weights(elph_ds,elph_tr_ds)
    
 use defs_basis
 use defs_elphon
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'd2c_weights'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(elph_type),intent(inout) :: elph_ds
 type(elph_tr_type),intent(inout),optional :: elph_tr_ds

!Local variables-------------------------------
 integer :: ii, jj, kk 
 integer :: ikpt, jkpt, kkpt
 integer :: iikpt, jjkpt, kkkpt
 integer :: icomp, jcomp
 integer :: iFSband, iband
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
 
!bxu, get wtk of coarse grid from fine grid
 elph_ds%k_phon%wtk = zero
 if (present(elph_tr_ds)) then
   elph_ds%k_phon%velocwtk = zero
   elph_ds%k_phon%vvelocwtk = zero
 end if

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
             elph_ds%k_phon%wtk(:,ikpt_phon,:)=elph_ds%k_phon%wtk(:,ikpt_phon,:)+&
&             elph_ds%k_fine%wtk(:,ikpt_fine,:)
             if (present(elph_tr_ds)) then
               do iFSband=1,elph_ds%ngkkband !FS bands
                 iband=iFSband+elph_ds%minFSband-1 ! full bands
                 do icomp = 1, 3
                   elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)=elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)+&
&                   elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)
                   do jcomp = 1, 3
                     elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:) = &
&                     elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:)+&
&                     elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)* &
&                     elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,jcomp,:)
                   end do
                 end do
               end do
             end if
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
             elph_ds%k_phon%wtk(:,ikpt_phon,:)=elph_ds%k_phon%wtk(:,ikpt_phon,:)+&
&             0.5_dp*elph_ds%k_fine%wtk(:,ikpt_fine,:)
             if (present(elph_tr_ds)) then
               do iFSband=1,elph_ds%ngkkband !FS bands
                 iband=iFSband+elph_ds%minFSband-1 ! full bands
                 do icomp = 1, 3
                   elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)=elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)+&
&                   0.5_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)
                   do jcomp = 1, 3
                     elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:) = &
&                     elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:)+&
&                     0.5_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)* &
&                     elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,jcomp,:)
                   end do
                 end do
               end do
             end if

             kkkpt = 1 + (kkpt-1)*nscale3 - nscale3/2
             if (kkkpt .le. 0) kkkpt = kkkpt + nkpt_fine3
             if (kkkpt .gt. nkpt_fine3) kkkpt = kkkpt - nkpt_fine3
             ikpt_fine = kkkpt + (jjkpt-1)*nkpt_fine3 + (iikpt-1)*nkpt_fine2*nkpt_fine3
             elph_ds%k_phon%wtk(:,ikpt_phon,:)=elph_ds%k_phon%wtk(:,ikpt_phon,:)+&
&             0.5_dp*elph_ds%k_fine%wtk(:,ikpt_fine,:)
             if (present(elph_tr_ds)) then
               do iFSband=1,elph_ds%ngkkband !FS bands
                 iband=iFSband+elph_ds%minFSband-1 ! full bands
                 do icomp = 1, 3
                   elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)=elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)+&
&                   0.5_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)
                   do jcomp = 1, 3
                     elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:) = &
&                     elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:)+&
&                     0.5_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)* &
&                     elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,jcomp,:)
                   end do
                 end do
               end do
             end if
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
             elph_ds%k_phon%wtk(:,ikpt_phon,:)=elph_ds%k_phon%wtk(:,ikpt_phon,:)+&
&             0.5_dp*elph_ds%k_fine%wtk(:,ikpt_fine,:)
             if (present(elph_tr_ds)) then
               do iFSband=1,elph_ds%ngkkband !FS bands
                 iband=iFSband+elph_ds%minFSband-1 ! full bands
                 do icomp = 1, 3
                   elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)=elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)+&
&                   0.5_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)
                   do jcomp = 1, 3
                     elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:) = &
&                     elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:)+&
&                     0.5_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)* &
&                     elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,jcomp,:)
                   end do
                 end do
               end do
             end if

             jjkpt = 1 + (jkpt-1)*nscale2 - nscale2/2
             if (jjkpt .le. 0) jjkpt = jjkpt + nkpt_fine2
             if (jjkpt .gt. nkpt_fine2) jjkpt = jjkpt - nkpt_fine2
             ikpt_fine = kkkpt + (jjkpt-1)*nkpt_fine3 + (iikpt-1)*nkpt_fine2*nkpt_fine3
             elph_ds%k_phon%wtk(:,ikpt_phon,:)=elph_ds%k_phon%wtk(:,ikpt_phon,:)+&
&             0.5_dp*elph_ds%k_fine%wtk(:,ikpt_fine,:)
             if (present(elph_tr_ds)) then
               do iFSband=1,elph_ds%ngkkband !FS bands
                 iband=iFSband+elph_ds%minFSband-1 ! full bands
                 do icomp = 1, 3
                   elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)=elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)+&
&                   0.5_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)
                   do jcomp = 1, 3
                     elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:) = &
&                     elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:)+&
&                     0.5_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)* &
&                     elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,jcomp,:)
                   end do
                 end do
               end do
             end if
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
             elph_ds%k_phon%wtk(:,ikpt_phon,:)=elph_ds%k_phon%wtk(:,ikpt_phon,:)+&
&             0.5_dp*elph_ds%k_fine%wtk(:,ikpt_fine,:)
             if (present(elph_tr_ds)) then
               do iFSband=1,elph_ds%ngkkband !FS bands
                 iband=iFSband+elph_ds%minFSband-1 ! full bands
                 do icomp = 1, 3
                   elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)=elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)+&
&                   0.5_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)
                   do jcomp = 1, 3
                     elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:) = &
&                     elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:)+&
&                     0.5_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)* &
&                     elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,jcomp,:)
                   end do
                 end do
               end do
             end if

             iikpt = 1 + (ikpt-1)*nscale1 - nscale1/2
             if (iikpt .le. 0) iikpt = iikpt + nkpt_fine1
             if (iikpt .gt. nkpt_fine1) iikpt = iikpt - nkpt_fine1
             ikpt_fine = kkkpt + (jjkpt-1)*nkpt_fine3 + (iikpt-1)*nkpt_fine2*nkpt_fine3
             elph_ds%k_phon%wtk(:,ikpt_phon,:)=elph_ds%k_phon%wtk(:,ikpt_phon,:)+&
&             0.5_dp*elph_ds%k_fine%wtk(:,ikpt_fine,:)
             if (present(elph_tr_ds)) then
               do iFSband=1,elph_ds%ngkkband !FS bands
                 iband=iFSband+elph_ds%minFSband-1 ! full bands
                 do icomp = 1, 3
                   elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)=elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)+&
&                   0.5_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)
                   do jcomp = 1, 3
                     elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:) = &
&                     elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:)+&
&                     0.5_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)* &
&                     elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,jcomp,:)
                   end do
                 end do
               end do
             end if
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
           elph_ds%k_phon%wtk(:,ikpt_phon,:)=elph_ds%k_phon%wtk(:,ikpt_phon,:)+&
&           0.25_dp*elph_ds%k_fine%wtk(:,ikpt_fine,:)
           if (present(elph_tr_ds)) then
             do iFSband=1,elph_ds%ngkkband !FS bands
               iband=iFSband+elph_ds%minFSband-1 ! full bands
               do icomp = 1, 3
                 elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)=elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)+&
&                 0.25_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)
                 do jcomp = 1, 3
                   elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:) = &
&                   elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:)+&
&                   0.25_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)* &
&                   elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,jcomp,:)
                 end do
               end do
             end do
           end if

           jjkpt = 1 + (jkpt-1)*nscale2 + nscale2/2
           kkkpt = 1 + (kkpt-1)*nscale3 - nscale3/2
           if (jjkpt .le. 0) jjkpt = jjkpt + nkpt_fine2
           if (jjkpt .gt. nkpt_fine2) jjkpt = jjkpt - nkpt_fine2
           if (kkkpt .le. 0) kkkpt = kkkpt + nkpt_fine3
           if (kkkpt .gt. nkpt_fine3) kkkpt = kkkpt - nkpt_fine3
           ikpt_fine = kkkpt + (jjkpt-1)*nkpt_fine3 + (iikpt-1)*nkpt_fine2*nkpt_fine3
           elph_ds%k_phon%wtk(:,ikpt_phon,:)=elph_ds%k_phon%wtk(:,ikpt_phon,:)+&
&           0.25_dp*elph_ds%k_fine%wtk(:,ikpt_fine,:)
           if (present(elph_tr_ds)) then
             do iFSband=1,elph_ds%ngkkband !FS bands
               iband=iFSband+elph_ds%minFSband-1 ! full bands
               do icomp = 1, 3
                 elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)=elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)+&
&                 0.25_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)
                 do jcomp = 1, 3
                   elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:) = &
&                   elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:)+&
&                   0.25_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)* &
&                   elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,jcomp,:)
                 end do
               end do
             end do
           end if

           jjkpt = 1 + (jkpt-1)*nscale2 - nscale2/2
           kkkpt = 1 + (kkpt-1)*nscale3 + nscale3/2
           if (jjkpt .le. 0) jjkpt = jjkpt + nkpt_fine2
           if (jjkpt .gt. nkpt_fine2) jjkpt = jjkpt - nkpt_fine2
           if (kkkpt .le. 0) kkkpt = kkkpt + nkpt_fine3
           if (kkkpt .gt. nkpt_fine3) kkkpt = kkkpt - nkpt_fine3
           ikpt_fine = kkkpt + (jjkpt-1)*nkpt_fine3 + (iikpt-1)*nkpt_fine2*nkpt_fine3
           elph_ds%k_phon%wtk(:,ikpt_phon,:)=elph_ds%k_phon%wtk(:,ikpt_phon,:)+&
&           0.25_dp*elph_ds%k_fine%wtk(:,ikpt_fine,:)
           if (present(elph_tr_ds)) then
             do iFSband=1,elph_ds%ngkkband !FS bands
               iband=iFSband+elph_ds%minFSband-1 ! full bands
               do icomp = 1, 3
                 elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)=elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)+&
&                 0.25_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)
                 do jcomp = 1, 3
                   elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:) = &
&                   elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:)+&
&                   0.25_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)* &
&                   elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,jcomp,:)
                 end do
               end do
             end do
           end if

           jjkpt = 1 + (jkpt-1)*nscale2 - nscale2/2
           kkkpt = 1 + (kkpt-1)*nscale3 - nscale3/2
           if (jjkpt .le. 0) jjkpt = jjkpt + nkpt_fine2
           if (jjkpt .gt. nkpt_fine2) jjkpt = jjkpt - nkpt_fine2
           if (kkkpt .le. 0) kkkpt = kkkpt + nkpt_fine3
           if (kkkpt .gt. nkpt_fine3) kkkpt = kkkpt - nkpt_fine3
           ikpt_fine = kkkpt + (jjkpt-1)*nkpt_fine3 + (iikpt-1)*nkpt_fine2*nkpt_fine3
           elph_ds%k_phon%wtk(:,ikpt_phon,:)=elph_ds%k_phon%wtk(:,ikpt_phon,:)+&
&           0.25_dp*elph_ds%k_fine%wtk(:,ikpt_fine,:)
           if (present(elph_tr_ds)) then
             do iFSband=1,elph_ds%ngkkband !FS bands
               iband=iFSband+elph_ds%minFSband-1 ! full bands
               do icomp = 1, 3
                 elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)=elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)+&
&                 0.25_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)
                 do jcomp = 1, 3
                   elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:) = &
&                   elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:)+&
&                   0.25_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)* &
&                   elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,jcomp,:)
                 end do
               end do
             end do
           end if
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
           elph_ds%k_phon%wtk(:,ikpt_phon,:)=elph_ds%k_phon%wtk(:,ikpt_phon,:)+&
&           0.25_dp*elph_ds%k_fine%wtk(:,ikpt_fine,:)
           if (present(elph_tr_ds)) then
             do iFSband=1,elph_ds%ngkkband !FS bands
               iband=iFSband+elph_ds%minFSband-1 ! full bands
               do icomp = 1, 3
                 elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)=elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)+&
&                 0.25_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)
                 do jcomp = 1, 3
                   elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:) = &
&                   elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:)+&
&                   0.25_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)* &
&                   elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,jcomp,:)
                 end do
               end do
             end do
           end if

           iikpt = 1 + (ikpt-1)*nscale1 + nscale1/2
           kkkpt = 1 + (kkpt-1)*nscale3 - nscale3/2
           if (iikpt .le. 0) iikpt = iikpt + nkpt_fine1
           if (iikpt .gt. nkpt_fine1) iikpt = iikpt - nkpt_fine1
           if (kkkpt .le. 0) kkkpt = kkkpt + nkpt_fine3
           if (kkkpt .gt. nkpt_fine3) kkkpt = kkkpt - nkpt_fine3
           ikpt_fine = kkkpt + (jjkpt-1)*nkpt_fine3 + (iikpt-1)*nkpt_fine2*nkpt_fine3
           elph_ds%k_phon%wtk(:,ikpt_phon,:)=elph_ds%k_phon%wtk(:,ikpt_phon,:)+&
&           0.25_dp*elph_ds%k_fine%wtk(:,ikpt_fine,:)
           if (present(elph_tr_ds)) then
             do iFSband=1,elph_ds%ngkkband !FS bands
               iband=iFSband+elph_ds%minFSband-1 ! full bands
               do icomp = 1, 3
                 elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)=elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)+&
&                 0.25_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)
                 do jcomp = 1, 3
                   elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:) = &
&                   elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:)+&
&                   0.25_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)* &
&                   elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,jcomp,:)
                 end do
               end do
             end do
           end if

           iikpt = 1 + (ikpt-1)*nscale1 - nscale1/2
           kkkpt = 1 + (kkpt-1)*nscale3 + nscale3/2
           if (iikpt .le. 0) iikpt = iikpt + nkpt_fine1
           if (iikpt .gt. nkpt_fine1) iikpt = iikpt - nkpt_fine1
           if (kkkpt .le. 0) kkkpt = kkkpt + nkpt_fine3
           if (kkkpt .gt. nkpt_fine3) kkkpt = kkkpt - nkpt_fine3
           ikpt_fine = kkkpt + (jjkpt-1)*nkpt_fine3 + (iikpt-1)*nkpt_fine2*nkpt_fine3
           elph_ds%k_phon%wtk(:,ikpt_phon,:)=elph_ds%k_phon%wtk(:,ikpt_phon,:)+&
&           0.25_dp*elph_ds%k_fine%wtk(:,ikpt_fine,:)
           if (present(elph_tr_ds)) then
             do iFSband=1,elph_ds%ngkkband !FS bands
               iband=iFSband+elph_ds%minFSband-1 ! full bands
               do icomp = 1, 3
                 elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)=elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)+&
&                 0.25_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)
                 do jcomp = 1, 3
                   elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:) = &
&                   elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:)+&
&                   0.25_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)* &
&                   elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,jcomp,:)
                 end do
               end do
             end do
           end if

           iikpt = 1 + (ikpt-1)*nscale1 - nscale1/2
           kkkpt = 1 + (kkpt-1)*nscale3 - nscale3/2
           if (iikpt .le. 0) iikpt = iikpt + nkpt_fine1
           if (iikpt .gt. nkpt_fine1) iikpt = iikpt - nkpt_fine1
           if (kkkpt .le. 0) kkkpt = kkkpt + nkpt_fine3
           if (kkkpt .gt. nkpt_fine3) kkkpt = kkkpt - nkpt_fine3
           ikpt_fine = kkkpt + (jjkpt-1)*nkpt_fine3 + (iikpt-1)*nkpt_fine2*nkpt_fine3
           elph_ds%k_phon%wtk(:,ikpt_phon,:)=elph_ds%k_phon%wtk(:,ikpt_phon,:)+&
&           0.25_dp*elph_ds%k_fine%wtk(:,ikpt_fine,:)
           if (present(elph_tr_ds)) then
             do iFSband=1,elph_ds%ngkkband !FS bands
               iband=iFSband+elph_ds%minFSband-1 ! full bands
               do icomp = 1, 3
                 elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)=elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)+&
&                 0.25_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)
                 do jcomp = 1, 3
                   elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:) = &
&                   elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:)+&
&                   0.25_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)* &
&                   elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,jcomp,:)
                 end do
               end do
             end do
           end if
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
           elph_ds%k_phon%wtk(:,ikpt_phon,:)=elph_ds%k_phon%wtk(:,ikpt_phon,:)+&
&           0.25_dp*elph_ds%k_fine%wtk(:,ikpt_fine,:)
           if (present(elph_tr_ds)) then
             do iFSband=1,elph_ds%ngkkband !FS bands
               iband=iFSband+elph_ds%minFSband-1 ! full bands
               do icomp = 1, 3
                 elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)=elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)+&
&                 0.25_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)
                 do jcomp = 1, 3
                   elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:) = &
&                   elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:)+&
&                   0.25_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)* &
&                   elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,jcomp,:)
                 end do
               end do
             end do
           end if

           jjkpt = 1 + (jkpt-1)*nscale2 + nscale2/2
           iikpt = 1 + (ikpt-1)*nscale1 - nscale1/2
           if (jjkpt .le. 0) jjkpt = jjkpt + nkpt_fine2
           if (jjkpt .gt. nkpt_fine2) jjkpt = jjkpt - nkpt_fine2
           if (iikpt .le. 0) iikpt = iikpt + nkpt_fine1
           if (iikpt .gt. nkpt_fine1) iikpt = iikpt - nkpt_fine1
           ikpt_fine = kkkpt + (jjkpt-1)*nkpt_fine3 + (iikpt-1)*nkpt_fine2*nkpt_fine3
           elph_ds%k_phon%wtk(:,ikpt_phon,:)=elph_ds%k_phon%wtk(:,ikpt_phon,:)+&
&           0.25_dp*elph_ds%k_fine%wtk(:,ikpt_fine,:)
           if (present(elph_tr_ds)) then
             do iFSband=1,elph_ds%ngkkband !FS bands
               iband=iFSband+elph_ds%minFSband-1 ! full bands
               do icomp = 1, 3
                 elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)=elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)+&
&                 0.25_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)
                 do jcomp = 1, 3
                   elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:) = &
&                   elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:)+&
&                   0.25_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)* &
&                   elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,jcomp,:)
                 end do
               end do
             end do
           end if

           jjkpt = 1 + (jkpt-1)*nscale2 - nscale2/2
           iikpt = 1 + (ikpt-1)*nscale1 + nscale1/2
           if (jjkpt .le. 0) jjkpt = jjkpt + nkpt_fine2
           if (jjkpt .gt. nkpt_fine2) jjkpt = jjkpt - nkpt_fine2
           if (iikpt .le. 0) iikpt = iikpt + nkpt_fine1
           if (iikpt .gt. nkpt_fine1) iikpt = iikpt - nkpt_fine1
           ikpt_fine = kkkpt + (jjkpt-1)*nkpt_fine3 + (iikpt-1)*nkpt_fine2*nkpt_fine3
           elph_ds%k_phon%wtk(:,ikpt_phon,:)=elph_ds%k_phon%wtk(:,ikpt_phon,:)+&
&           0.25_dp*elph_ds%k_fine%wtk(:,ikpt_fine,:)
           if (present(elph_tr_ds)) then
             do iFSband=1,elph_ds%ngkkband !FS bands
               iband=iFSband+elph_ds%minFSband-1 ! full bands
               do icomp = 1, 3
                 elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)=elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)+&
&                 0.25_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)
                 do jcomp = 1, 3
                   elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:) = &
&                   elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:)+&
&                   0.25_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)* &
&                   elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,jcomp,:)
                 end do
               end do
             end do
           end if

           jjkpt = 1 + (jkpt-1)*nscale2 - nscale2/2
           iikpt = 1 + (ikpt-1)*nscale1 - nscale1/2
           if (jjkpt .le. 0) jjkpt = jjkpt + nkpt_fine2
           if (jjkpt .gt. nkpt_fine2) jjkpt = jjkpt - nkpt_fine2
           if (iikpt .le. 0) iikpt = iikpt + nkpt_fine1
           if (iikpt .gt. nkpt_fine1) iikpt = iikpt - nkpt_fine1
           ikpt_fine = kkkpt + (jjkpt-1)*nkpt_fine3 + (iikpt-1)*nkpt_fine2*nkpt_fine3
           elph_ds%k_phon%wtk(:,ikpt_phon,:)=elph_ds%k_phon%wtk(:,ikpt_phon,:)+&
&           0.25_dp*elph_ds%k_fine%wtk(:,ikpt_fine,:)
           if (present(elph_tr_ds)) then
             do iFSband=1,elph_ds%ngkkband !FS bands
               iband=iFSband+elph_ds%minFSband-1 ! full bands
               do icomp = 1, 3
                 elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)=elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)+&
&                 0.25_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)
                 do jcomp = 1, 3
                   elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:) = &
&                   elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:)+&
&                   0.25_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)* &
&                   elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,jcomp,:)
                 end do
               end do
             end do
           end if
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
         elph_ds%k_phon%wtk(:,ikpt_phon,:)=elph_ds%k_phon%wtk(:,ikpt_phon,:)+&
&         0.125_dp*elph_ds%k_fine%wtk(:,ikpt_fine,:)
         if (present(elph_tr_ds)) then
           do iFSband=1,elph_ds%ngkkband !FS bands
             iband=iFSband+elph_ds%minFSband-1 ! full bands
             do icomp = 1, 3
               elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)=elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)+&
&               0.125_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)
               do jcomp = 1, 3
                 elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:) = &
&                 elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:)+&
&                 0.125_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)* &
&                 elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,jcomp,:)
               end do
             end do
           end do
         end if

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
         elph_ds%k_phon%wtk(:,ikpt_phon,:)=elph_ds%k_phon%wtk(:,ikpt_phon,:)+&
&         0.125_dp*elph_ds%k_fine%wtk(:,ikpt_fine,:)
         if (present(elph_tr_ds)) then
           do iFSband=1,elph_ds%ngkkband !FS bands
             iband=iFSband+elph_ds%minFSband-1 ! full bands
             do icomp = 1, 3
               elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)=elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)+&
&               0.125_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)
               do jcomp = 1, 3
                 elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:) = &
&                 elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:)+&
&                 0.125_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)* &
&                 elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,jcomp,:)
               end do
             end do
           end do
         end if

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
         elph_ds%k_phon%wtk(:,ikpt_phon,:)=elph_ds%k_phon%wtk(:,ikpt_phon,:)+&
&         0.125_dp*elph_ds%k_fine%wtk(:,ikpt_fine,:)
         if (present(elph_tr_ds)) then
           do iFSband=1,elph_ds%ngkkband !FS bands
             iband=iFSband+elph_ds%minFSband-1 ! full bands
             do icomp = 1, 3
               elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)=elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)+&
&               0.125_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)
               do jcomp = 1, 3
                 elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:) = &
&                 elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:)+&
&                 0.125_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)* &
&                 elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,jcomp,:)
               end do
             end do
           end do
         end if

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
         elph_ds%k_phon%wtk(:,ikpt_phon,:)=elph_ds%k_phon%wtk(:,ikpt_phon,:)+&
&         0.125_dp*elph_ds%k_fine%wtk(:,ikpt_fine,:)
         if (present(elph_tr_ds)) then
           do iFSband=1,elph_ds%ngkkband !FS bands
             iband=iFSband+elph_ds%minFSband-1 ! full bands
             do icomp = 1, 3
               elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)=elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)+&
&               0.125_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)
               do jcomp = 1, 3
                 elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:) = &
&                 elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:)+&
&                 0.125_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)* &
&                 elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,jcomp,:)
               end do
             end do
           end do
         end if

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
         elph_ds%k_phon%wtk(:,ikpt_phon,:)=elph_ds%k_phon%wtk(:,ikpt_phon,:)+&
&         0.125_dp*elph_ds%k_fine%wtk(:,ikpt_fine,:)
         if (present(elph_tr_ds)) then
           do iFSband=1,elph_ds%ngkkband !FS bands
             iband=iFSband+elph_ds%minFSband-1 ! full bands
             do icomp = 1, 3
               elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)=elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)+&
&               0.125_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)
               do jcomp = 1, 3
                 elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:) = &
&                 elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:)+&
&                 0.125_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)* &
&                 elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,jcomp,:)
               end do
             end do
           end do
         end if

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
         elph_ds%k_phon%wtk(:,ikpt_phon,:)=elph_ds%k_phon%wtk(:,ikpt_phon,:)+&
&         0.125_dp*elph_ds%k_fine%wtk(:,ikpt_fine,:)
         if (present(elph_tr_ds)) then
           do iFSband=1,elph_ds%ngkkband !FS bands
             iband=iFSband+elph_ds%minFSband-1 ! full bands
             do icomp = 1, 3
               elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)=elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)+&
&               0.125_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)
               do jcomp = 1, 3
                 elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:) = &
&                 elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:)+&
&                 0.125_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)* &
&                 elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,jcomp,:)
               end do
             end do
           end do
         end if

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
         elph_ds%k_phon%wtk(:,ikpt_phon,:)=elph_ds%k_phon%wtk(:,ikpt_phon,:)+&
&         0.125_dp*elph_ds%k_fine%wtk(:,ikpt_fine,:)
         if (present(elph_tr_ds)) then
           do iFSband=1,elph_ds%ngkkband !FS bands
             iband=iFSband+elph_ds%minFSband-1 ! full bands
             do icomp = 1, 3
               elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)=elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)+&
&               0.125_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)
               do jcomp = 1, 3
                 elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:) = &
&                 elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:)+&
&                 0.125_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)* &
&                 elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,jcomp,:)
               end do
             end do
           end do
         end if

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
         elph_ds%k_phon%wtk(:,ikpt_phon,:)=elph_ds%k_phon%wtk(:,ikpt_phon,:)+&
&         0.125_dp*elph_ds%k_fine%wtk(:,ikpt_fine,:)
         if (present(elph_tr_ds)) then
           do iFSband=1,elph_ds%ngkkband !FS bands
             iband=iFSband+elph_ds%minFSband-1 ! full bands
             do icomp = 1, 3
               elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)=elph_ds%k_phon%velocwtk(iFSband,ikpt_phon,icomp,:)+&
&               0.125_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)
               do jcomp = 1, 3
                 elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:) = &
&                 elph_ds%k_phon%vvelocwtk(iFSband,ikpt_phon,icomp,jcomp,:)+&
&                 0.125_dp*elph_ds%k_fine%wtk(iFSband,ikpt_fine,:)* &
&                 elph_tr_ds%el_veloc(ikpt_fine,iband,icomp,:)*elph_tr_ds%el_veloc(ikpt_fine,iband,jcomp,:)
               end do
             end do
           end do
         end if
       end if
     end do
   end do
 end do

!bxu, divide by nscale^3 to be consistent with the normalization of kpt_phon
 elph_ds%k_phon%wtk = elph_ds%k_phon%wtk/nscale1/nscale2/nscale3
 if (present(elph_tr_ds)) then
   elph_ds%k_phon%velocwtk = elph_ds%k_phon%velocwtk/nscale1/nscale2/nscale3
   elph_ds%k_phon%vvelocwtk = elph_ds%k_phon%vvelocwtk/nscale1/nscale2/nscale3
 end if

end subroutine d2c_weights
!!***

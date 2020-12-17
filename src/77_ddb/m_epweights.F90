!!****m* ABINIT/m_epweights
!! NAME
!!  m_epweights
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 2012-2020 ABINIT group (BXU, MVer)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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

module m_epweights

 use defs_basis
 use defs_elphon
 use m_abicore
 use m_errors
 use m_tetrahedron
 !use m_htetra
 use m_xmpi

 use m_symtk,           only : matr3inv

 implicit none

 private
!!***

 public :: d2c_weights
 public :: d2c_wtq!
 public :: ep_el_weights
 public :: ep_fs_weights
 public :: ep_ph_weights
!!***

contains
!!***

!!****f* ABINIT/d2c_weights
!! NAME
!!  d2c_weights
!!
!! FUNCTION
!!  This routine calculates the integration weights on a coarse k-grid
!!  using the integration weights from a denser k-gird. The weights of
!!  the extra k points that being shared are evenly distributed. It also
!!  condenses the velocity*wtk and velcity^2*wtk.
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
!!      m_elphon
!!
!! CHILDREN
!!      destroy_tetra,get_tetra_weight,init_tetra,matr3inv
!!
!! SOURCE

subroutine d2c_weights(elph_ds,elph_tr_ds)

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

!!****f* ABINIT/d2c_wtq
!! NAME
!!  d2c_wtq
!!
!! FUNCTION
!!  This routine calculates the integration weights on a coarse k-grid
!!  using the integration weights from a denser k-gird. The weights of
!!  the extra k points that being shared are evenly distributed.
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
!!      m_a2ftr,m_elphon
!!
!! CHILDREN
!!      destroy_tetra,get_tetra_weight,init_tetra,matr3inv
!!
!! SOURCE

subroutine d2c_wtq(elph_ds)

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

!!****f* ABINIT/ep_el_weights
!!
!! NAME
!! ep_el_weights
!!
!! FUNCTION
!! This routine calculates the Fermi Surface integration weights
!! for the electron phonon routines, by different methods
!!
!!    1) Gaussian smearing
!!    2) Tetrahedron method
!!    3) Window in bands for all k-points
!!    4) Fermi Dirac smearing, follows gaussian with a different smearing function
!!
!! INPUTS
!!   ep_b_min = minimal band to include in FS window integration
!!   ep_b_max = maximal band to include in FS window integration
!!   eigenGS = Ground State eigenvalues
!!   elphsmear = smearing width for Gaussian method
!!   fermie = Fermi level
!!   gprimd = Reciprocal lattice vectors (dimensionful)
!!   irredtoGS = mapping of elph k-points to ground state grid
!!   kptrlatt = k-point grid vectors (if divided by determinant of present matrix)
!!   max_occ = maximal occupancy for a band
!!   minFSband = minimal band included for Fermi Surface integration in Gaussian and Tetrahedron cases
!!   nFSband = number of bands in FS integration
!!   nsppol = number of spin polarizations
!!   telphint = option for FS integration:
!!      0 Tetrahedron method
!!      1 Gaussian smearing
!!      2 Window in bands for all k-points
!!      3 Fermi Dirac smearing
!!   k_obj%nkpt = number of FS k-points
!!   k_obj%kpt = FS k-points
!!   k_obj%full2irr = mapping of FS k-points from full grid to irred points
!!   k_obj%full2full = mapping of FS k-points in full grid under symops
!!
!! OUTPUT
!!
!! TODO
!!   weights should be recalculated on-the-fly! The present implementation is not flexible!
!!
!! PARENTS
!!      m_a2ftr,m_elphon
!!
!! CHILDREN
!!      destroy_tetra,get_tetra_weight,init_tetra,matr3inv
!!
!! SOURCE

subroutine ep_el_weights(ep_b_min, ep_b_max, eigenGS, elphsmear, enemin, enemax, nene, gprimd, &
&    irredtoGS, kptrlatt, max_occ, minFSband, nband, nFSband, nsppol, telphint, k_obj, tmp_wtk)

!Arguments ------------------------------------
!scalars
 type(elph_kgrid_type), intent(in) :: k_obj
 integer, intent(in) :: ep_b_min
 integer, intent(in) :: ep_b_max
 integer,intent(in) :: nband,nene
 real(dp), intent(in) :: elphsmear
 real(dp), intent(in) :: enemin,enemax
 real(dp), intent(in) :: gprimd(3,3)
 integer, intent(in) :: kptrlatt(3,3)
 real(dp), intent(in) :: max_occ
 integer, intent(in) :: minFSband
 integer, intent(in) :: nFSband
 integer, intent(in) :: nsppol
 integer, intent(in) :: telphint

! arrays
 real(dp), intent(in) :: eigenGS(nband,k_obj%nkptirr,nsppol)
 real(dp), intent(out) :: tmp_wtk(nFSband,k_obj%nkpt,nsppol,nene)
 integer, intent(in) :: irredtoGS(k_obj%nkptirr)

!Local variables-------------------------------
!scalars
 integer,parameter :: bcorr0=0
 integer :: ikpt, ikptgs, ib1, iband
 integer :: ierr, ie, isppol
 real(dp) :: deltaene, rcvol, fermie
 real(dp) :: smdeltaprefactor, smdeltafactor, xx

! arrays
 real(dp) :: rlatt(3,3), klatt(3,3)
 real(dp), allocatable :: tmp_eigen(:), tweight(:,:), dtweightde(:,:)
 character (len=500) :: message
 character (len=80) :: errstr
 type(t_tetrahedron) :: tetrahedra
 !type(htetra_t) :: tetrahedra

! *************************************************************************

 ! Initialize tmp_wtk with zeros
 tmp_wtk = zero

 !write(std_out,*) 'ep_el : nkpt ', k_obj%nkpt
!===================================
!Set up integration weights for FS
!===================================
 deltaene = (enemax-enemin)/dble(nene-1)

 if (telphint == 0) then

!  =========================
!  Tetrahedron integration
!  =========================

   rlatt(:,:) = kptrlatt(:,:)
   call matr3inv(rlatt,klatt)

   call init_tetra(k_obj%full2full(1,1,:), gprimd,klatt,k_obj%kpt, k_obj%nkpt,&
     tetrahedra, ierr, errstr, xmpi_comm_self)
   !call htetra_init(tetra, k_obj%full2full(1,1,:), gprimd, klatt, k_obj%kpt, k_obj%nkpt, &
   !                 k_obk%nkptirr, nkpt_ibz, ierr, errstr, xmpi_comm_self)

   ABI_CHECK(ierr==0,errstr)

   rcvol = abs (gprimd(1,1)*(gprimd(2,2)*gprimd(3,3)-gprimd(3,2)*gprimd(2,3)) &
&   -gprimd(2,1)*(gprimd(1,2)*gprimd(3,3)-gprimd(3,2)*gprimd(1,3)) &
&   +gprimd(3,1)*(gprimd(1,2)*gprimd(2,3)-gprimd(2,2)*gprimd(1,3)))

!  fix small window around fermie for tetrahedron weight calculation
   deltaene = (enemax-enemin)/dble(nene-1)

   ABI_ALLOCATE(tmp_eigen,(k_obj%nkpt))
   ABI_ALLOCATE(tweight,(k_obj%nkpt,nene))
   ABI_ALLOCATE(dtweightde,(k_obj%nkpt,nene))

   do iband = 1,nFSband
!    for each spin pol
     do isppol=1,nsppol
!    For this band get its contribution
       tmp_eigen(:) = zero
       do ikpt=1,k_obj%nkpt
         ikptgs = irredtoGS(k_obj%full2irr(1,ikpt))
         tmp_eigen(ikpt) = eigenGS(minFSband+iband-1,ikptgs,isppol)
       end do
!      calculate general integration weights at each irred kpoint
!      as in Blochl et al PRB 49 16223 [[cite:Bloechl1994a]]
       call get_tetra_weight(tmp_eigen,enemin,enemax,&
&       max_occ,nene,k_obj%nkpt,tetrahedra,bcorr0,&
&       tweight,dtweightde,xmpi_comm_self)

       tmp_wtk(iband,:,isppol,:) = dtweightde(:,:)*k_obj%nkpt
     end do
   end do
   ABI_DEALLOCATE(tmp_eigen)
   ABI_DEALLOCATE(tweight)
   ABI_DEALLOCATE(dtweightde)

   call destroy_tetra(tetrahedra)
   !call tetrahedra%free()

 else if (telphint == 1) then

!  ==============================================================
!  Gaussian or integration:
!  Each kpt contributes a gaussian of integrated weight 1
!  for each band. The gaussian being centered at the Fermi level.
!  ===============================================================

!  took out factor 1/k_obj%nkpt which intervenes only at integration time

!  MJV 18/5/2008 does smdeltaprefactor need to contain max_occ?

!  gaussian smdeltaprefactor = sqrt(piinv)/elphsmear/k_obj%nkpt
   smdeltaprefactor = max_occ*sqrt(piinv)/elphsmear
   smdeltafactor = one/elphsmear

!  SPPOL loop on isppol as well to get 2 sets of weights
   do isppol=1,nsppol
     fermie = enemin
     do ie = 1, nene
       fermie = fermie + deltaene
!      fine grid
       do ikpt=1, k_obj%nkpt
         ikptgs = irredtoGS(k_obj%full2irr(1,ikpt))
         do ib1=1,nFSband
           xx = smdeltafactor*(eigenGS(minFSband-1+ib1,ikptgs,isppol) - fermie)
           if (abs(xx) < 40._dp) then
             tmp_wtk(ib1,ikpt,isppol,ie) = exp(-xx*xx)*smdeltaprefactor
           end if
         end do
       end do
     end do
   end do


 else if (telphint == 2) then ! range of bands occupied

!  SPPOL eventually be able to specify bands for up and down separately
   fermie = enemin
   do ie = 1, nene
     fermie = fermie + deltaene
     do ikpt=1,k_obj%nkpt
       do ib1=ep_b_min, ep_b_max
!        for the moment both spin channels same
         tmp_wtk(ib1,ikpt,:,ie) = max_occ
       end do
     end do
   end do

   write(std_out,*) ' ep_el_weights : DOS is calculated from states in bands ',ep_b_min,' to ',ep_b_max

 else if (telphint == 3) then

!  ==============================================================
!  Fermi Dirac integration:
!  Each kpt contributes a Fermi Dirac smearing function of integrated weight 1
!  for each band. The function being centered at the Fermi level.
!  ===============================================================

!  took out factor 1/k_obj%nkpt which intervenes only at integration time

!  MJV 18/5/2008 does smdeltaprefactor need to contain max_occ?

!  gaussian smdeltaprefactor = sqrt(piinv)/elphsmear/k_obj%nkpt
   smdeltaprefactor = half*max_occ/elphsmear
   smdeltafactor = one/elphsmear

!  SPPOL loop on isppol as well to get 2 sets of weights
   do isppol=1,nsppol
     fermie = enemin
     do ie = 1, nene
       fermie = fermie + deltaene
!      fine grid
       do ikpt=1, k_obj%nkpt
         ikptgs = irredtoGS(k_obj%full2irr(1,ikpt))
         do ib1=1,nFSband
           xx = smdeltafactor*(eigenGS(minFSband-1+ib1,ikptgs,isppol) - fermie)
           tmp_wtk(ib1,ikpt,isppol,ie) = smdeltaprefactor / (one + cosh(xx))
         end do
       end do
     end do
   end do


 else
   write (message,'(a,i0)')" telphint should be between 0 and 3, found: ",telphint
   MSG_BUG(message)
 end if ! if telphint

end subroutine ep_el_weights
!!***

!!****f* ABINIT/ep_fs_weights
!!
!! NAME
!! ep_fs_weights
!!
!! FUNCTION
!! This routine calculates the Fermi Surface integration weights
!!  for the electron phonon routines, by different methods
!!    1) Gaussian smearing
!!    2) Tetrahedron method
!!    3) Window in bands for all k-points
!!    4) Fermi Dirac smearing, follows gaussian with a different smearing function
!!
!! INPUTS
!!   ep_b_min = minimal band to include in FS window integration
!!   ep_b_max = maximal band to include in FS window integration
!!   eigenGS = Ground State eigenvalues
!!   elphsmear = smearing width for Gaussian method
!!   fermie = Fermi level
!!   gprimd = Reciprocal lattice vectors (dimensionful)
!!   irredtoGS = mapping of elph k-points to ground state grid
!!   kptrlatt = k-point grid vectors (if divided by determinant of present matrix)
!!   max_occ = maximal occupancy for a band
!!   minFSband = minimal band included for Fermi Surface integration in Gaussian and Tetrahedron cases
!!   nFSband = number of bands in FS integration
!!   nsppol = number of spin polarizations
!!   telphint = option for FS integration:
!!      0 Tetrahedron method
!!      1 Gaussian smearing
!!      2 Window in bands for all k-points
!!      3 Fermi Dirac smearing
!!   k_obj%nkpt = number of FS k-points
!!   k_obj%kpt = FS k-points
!!   k_obj%full2irr = mapping of FS k-points from full grid to irred points
!!   k_obj%full2full = mapping of FS k-points in full grid under symops
!!
!! OUTPUT
!!   k_obj%wtk = integration weights
!!
!! TODO
!!   weights should be recalculated on-the-fly! The present implementation is not flexible!
!!
!! PARENTS
!!      m_elphon
!!
!! CHILDREN
!!      destroy_tetra,get_tetra_weight,init_tetra,matr3inv
!!
!! SOURCE

subroutine ep_fs_weights(ep_b_min, ep_b_max, eigenGS, elphsmear, fermie, gprimd, &
&    irredtoGS, kptrlatt, max_occ, minFSband, nband, nFSband, nsppol, telphint, k_obj)

!Arguments ------------------------------------
!scalars
 type(elph_kgrid_type), intent(inout) :: k_obj
 integer, intent(in) :: ep_b_min
 integer, intent(in) :: ep_b_max
 integer,intent(in) :: nband
 real(dp), intent(in) :: elphsmear
 real(dp), intent(in) :: fermie
 real(dp), intent(in) :: gprimd(3,3)
 integer, intent(in) :: kptrlatt(3,3)
 real(dp), intent(in) :: max_occ
 integer, intent(in) :: minFSband
 integer, intent(in) :: nFSband
 integer, intent(in) :: nsppol
 integer, intent(in) :: telphint

! arrays
 real(dp), intent(in) :: eigenGS(nband,k_obj%nkptirr,nsppol)
 integer, intent(in) :: irredtoGS(k_obj%nkptirr)

!Local variables-------------------------------
!scalars
 integer,parameter :: bcorr0=0
 integer :: ikpt, ikptgs, ib1, isppol, iband
 integer :: nene, ifermi
 integer :: ierr

 real(dp) :: enemin, enemax, deltaene, rcvol
 real(dp) :: smdeltaprefactor, smdeltafactor, xx

! arrays
 real(dp) :: rlatt(3,3), klatt(3,3)
 real(dp), allocatable :: tmp_eigen(:), tweight(:,:), dtweightde(:,:)

 character (len=500) :: message
 character (len=80) :: errstr

 type(t_tetrahedron) :: tetrahedra

! *************************************************************************

 write(std_out,*) 'ep_fs : nkpt ', k_obj%nkpt
 write(message, '(a)' ) '- ep_fs_weights  1  = '
 call wrtout(std_out,message,'PERS')

!===================================
!Set up integration weights for FS
!===================================

 if (telphint == 0) then

!  =========================
!  Tetrahedron integration
!  =========================

   rlatt(:,:) = kptrlatt(:,:)
   call matr3inv(rlatt,klatt)

   call init_tetra(k_obj%full2full(1,1,:), gprimd,klatt,k_obj%kpt, k_obj%nkpt,&
&   tetrahedra, ierr, errstr, xmpi_comm_self)
   ABI_CHECK(ierr==0,errstr)

   rcvol = abs (gprimd(1,1)*(gprimd(2,2)*gprimd(3,3)-gprimd(3,2)*gprimd(2,3)) &
&   -gprimd(2,1)*(gprimd(1,2)*gprimd(3,3)-gprimd(3,2)*gprimd(1,3)) &
&   +gprimd(3,1)*(gprimd(1,2)*gprimd(2,3)-gprimd(2,2)*gprimd(1,3)))

!  just do weights at FS
   nene = 100

!  fix small window around fermie for tetrahedron weight calculation
   deltaene = 2*elphsmear/dble(nene-1)
   ifermi = int(nene/2)
   enemin = fermie - dble(ifermi-1)*deltaene
   enemax = enemin + dble(nene-1)*deltaene

   ABI_ALLOCATE(tmp_eigen,(k_obj%nkpt))
   ABI_ALLOCATE(tweight,(k_obj%nkpt,nene))
   ABI_ALLOCATE(dtweightde,(k_obj%nkpt,nene))

   do iband = 1,nFSband
!    for each spin pol
     do isppol=1,nsppol
!      For this band get its contribution
       tmp_eigen(:) = zero
       do ikpt=1,k_obj%nkpt
         ikptgs = irredtoGS(k_obj%full2irr(1,ikpt))
         tmp_eigen(ikpt) = eigenGS(minFSband+iband-1,ikptgs,isppol)
       end do
!      calculate general integration weights at each irred kpoint
!      as in Blochl et al PRB 49 16223 [[cite:Bloechl1994a]]
       call get_tetra_weight(tmp_eigen,enemin,enemax,&
&       max_occ,nene,k_obj%nkpt,tetrahedra,bcorr0,&
&       tweight,dtweightde,xmpi_comm_self)

       k_obj%wtk(iband,:,isppol) = dtweightde(:,ifermi)*k_obj%nkpt
     end do

   end do
   ABI_DEALLOCATE(tmp_eigen)
   ABI_DEALLOCATE(tweight)
   ABI_DEALLOCATE(dtweightde)

   call destroy_tetra(tetrahedra)

 else if (telphint == 1) then

!  ==============================================================
!  Gaussian or integration:
!  Each kpt contributes a gaussian of integrated weight 1
!  for each band. The gaussian being centered at the Fermi level.
!  ===============================================================

!  took out factor 1/k_obj%nkpt which intervenes only at integration time

!  MJV 18/5/2008 does smdeltaprefactor need to contain max_occ?

!  gaussian smdeltaprefactor = sqrt(piinv)/elphsmear/k_obj%nkpt
   smdeltaprefactor = max_occ*sqrt(piinv)/elphsmear
   smdeltafactor = one/elphsmear

   k_obj%wtk = zero
!  SPPOL loop on isppol as well to get 2 sets of weights
   do isppol=1,nsppol
!    fine grid
     do ikpt=1, k_obj%nkpt
       ikptgs = irredtoGS(k_obj%full2irr(1,ikpt))
       do ib1=1,nFSband
         xx = smdeltafactor*(eigenGS(minFSband-1+ib1,ikptgs,isppol) - fermie)
         if (abs(xx) < 40._dp) then
           k_obj%wtk(ib1,ikpt,isppol) = exp(-xx*xx)*smdeltaprefactor
         end if
       end do
     end do
   end do


 else if (telphint == 2) then ! range of bands occupied

!  SPPOL eventually be able to specify bands for up and down separately
   k_obj%wtk = zero
   do ikpt=1,k_obj%nkpt
     do ib1=ep_b_min, ep_b_max
!      for the moment both spin channels same
       k_obj%wtk(ib1,ikpt,:) = max_occ
     end do
   end do

   write(std_out,*) ' ep_fs_weights : DOS is calculated from states in bands ',ep_b_min,' to ',ep_b_max

 else if (telphint == 3) then

!  ==============================================================
!  Fermi Dirac integration:
!  Each kpt contributes a Fermi Dirac smearing function of integrated weight 1
!  for each band. The function being centered at the Fermi level.
!  ===============================================================

!  took out factor 1/k_obj%nkpt which intervenes only at integration time

!  MJV 18/5/2008 does smdeltaprefactor need to contain max_occ?

!  gaussian smdeltaprefactor = sqrt(piinv)/elphsmear/k_obj%nkpt
   smdeltaprefactor = half*max_occ/elphsmear
   smdeltafactor = one/elphsmear

   k_obj%wtk = zero
!  SPPOL loop on isppol as well to get 2 sets of weights
   do isppol=1,nsppol
!    fine grid
     do ikpt=1, k_obj%nkpt
       ikptgs = irredtoGS(k_obj%full2irr(1,ikpt))
       do ib1=1,nFSband
         xx = smdeltafactor*(eigenGS(minFSband-1+ib1,ikptgs,isppol) - fermie)
         k_obj%wtk(ib1,ikpt,isppol) = smdeltaprefactor / (one + cosh(xx))
       end do
     end do
   end do


 else
   write (message,'(a,i0)')" telphint should be between 0 and 3, found: ",telphint
   MSG_BUG(message)
 end if ! if telphint

end subroutine ep_fs_weights
!!***

!!****f* ABINIT/ep_ph_weights
!!
!! NAME
!! ep_ph_weights
!!
!! FUNCTION
!! This routine calculates the phonon integration weights
!!  for the electron phonon routines, by different methods
!!    1) Gaussian smearing
!!    0) Tetrahedron method
!!
!! INPUTS
!!   phfrq = phonon energies
!!   elphsmear = smearing width for Gaussian method
!!   omega = input phonon energy
!!   gprimd = Reciprocal lattice vectors (dimensionful)
!!   kptrlatt = k-point grid vectors (if divided by determinant of present matrix)
!!   telphint = option for FS integration:
!!      0 Tetrahedron method
!!      1 Gaussian smearing
!!   k_obj%nkpt = number of FS k-points
!!   k_obj%kpt = FS k-points
!!   k_obj%full2full = mapping of FS k-points in full grid under symops
!!
!! OUTPUT
!!   tmp_wtq = integration weights
!!
!! TODO
!!   weights should be recalculated on-the-fly! The present implementation is not flexible!
!!
!! PARENTS
!!      m_a2ftr,m_elphon
!!
!! CHILDREN
!!      destroy_tetra,get_tetra_weight,init_tetra,matr3inv
!!
!! SOURCE

subroutine ep_ph_weights(phfrq,elphsmear,omega_min,omega_max,nomega,gprimd,kptrlatt,nbranch,telphint,k_obj,tmp_wtq)

!Arguments ------------------------------------
!scalars
 type(elph_kgrid_type), intent(inout) :: k_obj
 integer,intent(in) :: nbranch
 real(dp), intent(in) :: elphsmear
 real(dp), intent(in) :: omega_min,omega_max
 real(dp), intent(in) :: gprimd(3,3)
 integer, intent(in) :: kptrlatt(3,3)
 integer, intent(in) :: nomega
 integer, intent(in) :: telphint

! arrays
 real(dp), intent(in) :: phfrq(nbranch,k_obj%nkpt)
 real(dp), intent(out) :: tmp_wtq(nbranch,k_obj%nkpt,nomega)

!Local variables-------------------------------
!scalars
 integer,parameter :: bcorr0=0
 integer :: ikpt, ib1, ibranch
 integer :: ierr, iomega
 real(dp) :: rcvol, max_occ
 real(dp) :: smdeltaprefactor, smdeltafactor, xx, gaussmaxarg
 real(dp) :: domega,omega

! arrays
 real(dp) :: rlatt(3,3), klatt(3,3)
 real(dp), allocatable :: tweight(:,:), dtweightde(:,:)
 character (len=80) :: errstr
 type(t_tetrahedron) :: tetrahedra

! *************************************************************************

 !write(std_out,*) 'ep_ph : nqpt ', k_obj%nkpt
!===================================
!Set up integration weights for FS
!===================================
 max_occ = one
 gaussmaxarg = sqrt(-log(1.d-100))
 domega = (omega_max - omega_min)/(nomega - 1)

 if (telphint == 0) then

!  =========================
!  Tetrahedron integration
!  =========================

   rlatt(:,:) = kptrlatt(:,:)
   call matr3inv(rlatt,klatt)

   call init_tetra(k_obj%full2full(1,1,:), gprimd,klatt,k_obj%kpt, k_obj%nkpt,&
&   tetrahedra, ierr, errstr, xmpi_comm_self)
   ABI_CHECK(ierr==0,errstr)

   rcvol = abs (gprimd(1,1)*(gprimd(2,2)*gprimd(3,3)-gprimd(3,2)*gprimd(2,3)) &
&   -gprimd(2,1)*(gprimd(1,2)*gprimd(3,3)-gprimd(3,2)*gprimd(1,3)) &
&   +gprimd(3,1)*(gprimd(1,2)*gprimd(2,3)-gprimd(2,2)*gprimd(1,3)))

!  do all the omega points for tetrahedron weight calculation

   ABI_ALLOCATE(tweight,(k_obj%nkpt,nomega))
   ABI_ALLOCATE(dtweightde,(k_obj%nkpt,nomega))

   do ibranch = 1,nbranch
     call get_tetra_weight(phfrq(ibranch,:),omega_min,omega_max,&
&     max_occ,nomega,k_obj%nkpt,tetrahedra,bcorr0,&
&     tweight,dtweightde,xmpi_comm_self)

     tmp_wtq(ibranch,:,:) = dtweightde(:,:)*k_obj%nkpt
   end do
   ABI_DEALLOCATE(tweight)
   ABI_DEALLOCATE(dtweightde)

   call destroy_tetra(tetrahedra)

 else if (telphint == 1) then

!  ==============================================================
!  Gaussian or integration:
!  Each kpt contributes a gaussian of integrated weight 1
!  for each branch. The gaussian being centered at the input energy
!  ===============================================================

!  took out factor 1/k_obj%nkpt which intervenes only at integration time

!  gaussian smdeltaprefactor = sqrt(piinv)/elphsmear/k_obj%nkpt
   smdeltaprefactor = max_occ*sqrt(piinv)/elphsmear
   smdeltafactor = one/elphsmear

   tmp_wtq = zero
   omega = omega_min
   do iomega = 1, nomega
     omega = omega + domega
     do ikpt=1, k_obj%nkpt
       do ib1=1,nbranch
         xx = smdeltafactor*(phfrq(ib1,ikpt)-omega)
         if (abs(xx) < gaussmaxarg) then
           tmp_wtq(ib1,ikpt,iomega) = exp(-xx*xx)*smdeltaprefactor
         end if
       end do
     end do
   end do
 end if ! if telphint

end subroutine ep_ph_weights
!!***

end module m_epweights
!!***

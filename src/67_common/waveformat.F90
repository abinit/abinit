!{\src2tex{textfont=tt}}
!!****f* ABINIT/waveformat
!! NAME
!! waveformat
!!
!!
!! FUNCTION
!! This routine is to find the matched pairs of plane waves between
!! two neighbouring k points and load a new pw coefficients array cg_new
!! Was written first by Na Sai (thanks), but unfortunately without
!! any comment ...
!!
!! COPYRIGHT
!! Copyright (C) 2000-2016 ABINIT group (NSAI,XG,MKV)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)= input planewave coefficients, in case mkmem/=0
!!  cg_disk(2,mpw*nspinor*mband,2)= input planewave coefficients, in case mkmem==0
!!  cg_index(mband,nkpt_,nsppol)=index of wavefunction iband,ikpt,isppol in the array cg.
!!  dk(3)= step taken to the next k mesh point along the kberry direction (see also isgn)
!!  ii=(to be documented)
!!  ikpt=index of the first k-point in the reduced Brillouin zone
!!  ikpt_=index of the first k-point in the full Brillouin zone
!!  isgn=1 if dk(3) is connecting the k-points (ikpt_ and jkpt)
!!      =-1 if -dk(3) is connecting the k-points
!!  isppol=1 if spin-up, =2 if spin-down
!!  jj=(to be documented)
!!  jkpt=index of the second k-point in the reduced Brillouin zone
!!  jkpt_=index of the second k-point in the full Brillouin zone
!!  kg_kpt(:,:,:)= unpacked reduced planewave coordinates with subscript of
!!          planewave and k point
!!  kpt(3,nkpt)=reduced coordinates of k-point grid that samples the whole BZ
!!  kg_jl(3,mpw,2)=(to be documented)
!!  maxband/minband= control the minimum and maximum band calculated in the
!!           overlap matrix
!!  mband=maximum number of bands (dimension of several cg* arrays)
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mcg_disk=size of wave-functions array (cg_disk) =mpw*nspinor*mband
!!  mkmem= if 0, the wavefunctions are input in cg_disk, otherwise in cg
!!  mpw=maximum number of planewaves (dimension of several cg* arrays)
!!  nkpt=number of k points (full Brillouin zone !?!)
!!  nkpt_=number of k points (reduced Brillouin zone !?!)
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  shift_g_2(nkpt,nkpt)=non-zero if a G vector along kberry is needed to connect k points
!!  tr(2)=variable that changes k to -k
!!                              G to -G
!!                     $c_g$ to $c_g^*$ when time-reversal symetrie is used
!!
!! OUTPUT
!!  cg_new(2,mpw,maxband)=planewave coefficients transferred onto the
!!   set of planewaves at k
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      uderiv
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine waveformat(cg,cg_disk,cg_index,cg_new,dk,ii,ikpt,&
& ikpt_,isgn,isppol,jj,jkpt,jkpt_,kg_kpt,kpt,kg_jl,maxband,mband,mcg,mcg_disk,&
& minband,mkmem,mpw,nkpt,nkpt_,npwarr,nsppol,nspinor,shift_g_2,tr)

 use m_profiling_abi

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'waveformat'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ii,ikpt,ikpt_,isgn,isppol,jj,jkpt,jkpt_,maxband,mband,mcg,mcg_disk
 integer,intent(in) :: minband,mkmem,mpw,nkpt,nkpt_,nspinor,nsppol
!arrays
 integer,intent(in) :: cg_index(mband,nkpt_,nsppol),kg_jl(3,mpw,2)
 integer,intent(in) :: kg_kpt(3,mpw*nspinor,nkpt_),npwarr(nkpt_)
 real(dp),intent(in) :: cg(2,mcg)
 real(dp),intent(in) :: cg_disk(2,mcg_disk,2),dk(3),kpt(3,nkpt),tr(2)
 real(dp),intent(out) :: cg_new(2,mpw,maxband)
 logical,intent(in) :: shift_g_2(nkpt,nkpt)

!Local variables -------------------------
!scalars
 integer :: cg_index_iband,iband,ipw,jpw,nomatch,npw_k
 logical :: found_match
!arrays
 integer :: dg(3)

! ***********************************************************************

 npw_k=npwarr(ikpt)


 nomatch=0

!If there is no shift of G-vector between ikpt_ and jkpt_
 if(shift_g_2(ikpt_,jkpt_) .eqv. .false.) then

!  DEBUG
!  write(111,*)'pair', ikpt_,jkpt_,'noshift'
!  ENDDEBUG

!  If the original wavefunction is contained in cg_disk
   if(mkmem==0) then

     do ipw=1,npw_k

       found_match = .false.

       do jpw=1,npwarr(jkpt)
         if (sum(abs(tr(ii)*kg_jl(:,ipw,ii)-tr(jj)*kg_jl(:,jpw,jj)))<3*tol8)then
           do iband=minband, maxband
             cg_index_iband=(iband-1)*npwarr(jkpt)
             cg_new(1:2,ipw,iband)=cg_disk(1:2,jpw+cg_index_iband,jj)
           end do
           found_match = .true.
           exit
         end if
       end do

       if (found_match .eqv. .false.) then
         do iband=minband,maxband
           cg_new(1:2,ipw,iband)=zero
         end do
         nomatch = nomatch + 1
       end if

     end do

!    Here, the wavefunctions are contained in cg
   else

     do ipw=1,npw_k

       found_match = .false.

       do jpw=1,npwarr(jkpt)
         if (sum(abs(tr(ii)*kg_kpt(:,ipw,ikpt)-tr(jj)*kg_kpt(:,jpw,jkpt)))<3*tol8)then
           do iband=minband, maxband
             cg_index_iband=cg_index(iband,jkpt,isppol)
             cg_new(1:2,ipw,iband)=cg(1:2,jpw+cg_index_iband)
           end do
           found_match = .true.
           exit
         end if
       end do

       if (found_match .eqv. .false.) then
         do iband=minband,maxband
           cg_new(1:2,ipw,iband)=(0.0_dp,0.0_dp)
         end do
         nomatch = nomatch + 1
       end if
     end do

   end if

!  DEBUG
!  write(111,*) 'normal pair nomatch=',nomatch
!  ENDDEBUG

!  If there is a G-vector shift between ikpt_ and jkpt_
 else

!  DEBUG
!  write(111,*) 'pair',ikpt_,jkpt_,' need shift'
!  ENDDEBUG

   dg(:) = -1*nint(tr(jj)*kpt(:,jkpt)-tr(ii)*kpt(:,ikpt)+isgn*dk(:))

!  If the original wavefunction is contained in cg_disk
   if(mkmem==0) then

     do ipw=1,npw_k

       found_match = .false.

       do jpw=1,npwarr(jkpt)
         if (sum(abs(tr(ii)*kg_jl(:,ipw,ii)-(tr(jj)*kg_jl(:,jpw,jj)-&
&         dg(:))))<3*tol8)then

           do iband=minband, maxband
             cg_index_iband=(iband-1)*npwarr(jkpt)
             cg_new(1:2,ipw,iband)=cg_disk(1:2,jpw+cg_index_iband,jj)
           end do
           found_match = .true.
           exit
         end if
       end do

       if (found_match .eqv. .false.) then
         do iband=minband,maxband
           cg_new(1:2,ipw,iband)=(0.0_dp,0.0_dp)
         end do
         nomatch = nomatch + 1
       end if
     end do

!    Here, the wavefunctions are contained in cg
   else

     do ipw=1,npw_k

       found_match = .false.

       do jpw=1,npwarr(jkpt)
         if (sum(abs(tr(ii)*kg_kpt(:,ipw,ikpt)-(tr(jj)*kg_kpt(:,jpw,jkpt)-&
&         dg(:))))<3*tol8)then
           do iband=minband, maxband
             cg_index_iband=cg_index(iband,jkpt,isppol)
             cg_new(1:2,ipw,iband)=cg(1:2,jpw+cg_index_iband)
           end do
           found_match = .true.
           exit
         end if
       end do

       if (found_match .eqv. .false.) then
         do iband=minband,maxband
           cg_new(1:2,ipw,iband)=zero
         end do
         nomatch = nomatch + 1
       end if
     end do

   end if

!  DEBUG
!  write(111,*) 'special pair nomatch=',nomatch
!  ENDDEBUG

 end if

end subroutine waveformat
!!***

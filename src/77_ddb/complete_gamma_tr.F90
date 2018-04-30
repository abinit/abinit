!{\src2tex{textfont=tt}}
!!****f* ABINIT/complete_gamma_tr
!!
!! NAME
!! complete_gamma_tr
!!
!! FUNCTION
!! Use the set of special q points calculated by the Monkhorst & Pack Technique.
!! Check if all the informations for the q points are present in
!! the input gamma transport matrices.
!! Generate the gamma transport matrices (already summed over the FS) of the set of q points which
!! samples homogeneously the entire Brillouin zone.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2018 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! crystal<crystal_t>=data type gathering info on the crystalline structure.
!! ep_scalprod= flag for scalar product of gkk with phonon displacement vectors
!! nbranch=number of phonon branches = 3*natom
!! nqptirred=nqpt irred BZ
!! nqpt_full=nqpt full BZ
!! nsppol=number of spins
!! qirredtofull= mapping irred to full qpoints
!! qpttoqpt = qpoint index mapping under symops
!!
!! OUTPUT
!! gamma_qpt_tr = in/out: set of gamma matrix elements completed and symmetrized 
!!    gamma_qpt_tr(2,9,nbranch*nbranch,nsppol,nqpt_full)
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      dgemm
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine complete_gamma_tr(crystal,ep_scalprod,nbranch,nqptirred,nqpt_full,nsppol,gamma_qpt_tr,qirredtofull,qpttoqpt)

 use defs_basis
 use defs_elphon
 use m_profiling_abi
 use m_linalg_interfaces
 use m_errors

 use m_crystal,   only : crystal_t

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'complete_gamma_tr'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: nbranch,nqptirred,nqpt_full,nsppol, ep_scalprod
 type(crystal_t),intent(in) :: crystal
!arrays
 integer,intent(in) :: qpttoqpt(2,crystal%nsym,nqpt_full)
 integer,intent(in) :: qirredtofull(nqptirred)
 real(dp), intent(inout) :: gamma_qpt_tr(2,9,nbranch*nbranch,nsppol,nqpt_full)

!Local variables-------------------------------
!scalars
 integer :: ieqqpt,ii,iqpt,isppol,isym
 integer :: itim,jj,kk,ll,neqqpt
 integer :: iatom,ancestor_iatom
 integer :: iqpt_fullbz,imode, itensor
 real(dp),parameter :: tol=2.d-8
!arrays
 integer :: symrel(3,3,crystal%nsym),symrec(3,3,crystal%nsym)
 integer :: symmetrized_qpt(nqpt_full)
 integer :: gkk_flag(nbranch,nbranch,nsppol,nqpt_full)
 real(dp) :: gprimd(3,3),rprimd(3,3)
 real(dp) :: ss(3,3), sscart(3,3)
 real(dp) :: tmp_mat(2,nbranch,nbranch)
 real(dp) :: tmp_mat2(2,nbranch,nbranch)
 real(dp) :: tmp_tensor(2,3,3)
 real(dp) :: tmp_tensor2(2,3,3)
 real(dp) :: ss_allatoms(nbranch,nbranch)
 real(dp) :: c_one(2), c_zero(2)
 real(dp),allocatable :: gkk_qpt_new(:,:,:,:),gkk_qpt_tmp(:,:,:,:)

! *********************************************************************

 c_one = (/one,zero/)
 c_zero = (/zero,zero/)

 gprimd = crystal%gprimd
 rprimd = crystal%rprimd

 symrec =  crystal%symrec
 symrel =  crystal%symrel

!Generation of the gkk matrices relative to the q points
!of the set which samples the entire Brillouin zone

!set up flags for gamma_qpt matrices we have
 gkk_flag = -1
 do iqpt=1,nqptirred
   iqpt_fullbz = qirredtofull(iqpt)
   gkk_flag(:,:,:,iqpt_fullbz) = 1
 end do

 symmetrized_qpt(:) = -1
! isppol=1

 ABI_ALLOCATE(gkk_qpt_new,(2,9,nbranch*nbranch, nsppol))
 ABI_ALLOCATE(gkk_qpt_tmp,(2,9,nbranch*nbranch, nsppol))

 do iqpt=1,nqpt_full

!  Already symmetrized?
   if (symmetrized_qpt(iqpt) == 1) cycle

   gkk_qpt_new(:,:,:,:) = zero

!  loop over qpoints equivalent to iqpt
   neqqpt=0
!  do not use time reversal symmetry to complete the qpoints:
!  do not know what happens to the gamma matrices

   do itim=1,2
     do isym=1,crystal%nsym
!      ieqqpt is sent onto iqpt by itim/isym
       ieqqpt = qpttoqpt(itim,isym,iqpt)

       if (gkk_flag(1,1,1,ieqqpt) == -1) cycle
!      if we have information on this qpt
!      iqpt is equivalent to ieqqpt: get it from file or memory
       gkk_qpt_tmp(:,:,:,:) = gamma_qpt_tr(:,:,:,:,ieqqpt)

       neqqpt=neqqpt+1

!      
!      MJV note 02/2010:
!      the correspondence of symrel and symrec in the different cases, symmetrizing there
!      and back, has been fixed in the cases with and without scalprod (ie cartesian
!      and reduced real space coordinates) with respect to a calculation with no symmetries
!      I believe everything is settled, but still do not know why the 2 versions of the ss
!      matrices here use different rel/rec, instead of just being multiplied by the rprim gprim...
!      
       do ii=1,3
         do jj=1,3
           sscart(ii,jj)=0.0_dp
           do kk=1,3
             do ll=1,3
               sscart(ii,jj)=sscart(ii,jj)+rprimd(ii,kk)*symrel(kk,ll,isym)*gprimd(ll,jj)
!              sscart(ii,jj)=sscart(ii,jj)+rprimd(ii,kk)*symrel(kk,ll,isym)*gprimd(ll,jj)
             end do
           end do
         end do
       end do
       if (ep_scalprod==1) then
         do ii=1,3
           do jj=1,3
             ss(ii,jj)=0.0_dp
             do kk=1,3
               do ll=1,3
                 ss(ii,jj)=ss(ii,jj)+rprimd(ii,kk)*symrel(kk,ll,isym)*gprimd(ll,jj)
               end do
             end do
           end do
         end do
       else
         do ii=1,3
           do jj=1,3
             ss(ii,jj) = symrec(ii,jj,isym)
           end do
         end do
       end if

       ss_allatoms(:,:) = zero
       do iatom=1,crystal%natom
         ancestor_iatom = crystal%indsym(4,isym,iatom)
         ss_allatoms((ancestor_iatom-1)*3+1:(ancestor_iatom-1)*3+3,&
&         (iatom-1)*3+1:         (iatom-1)*3+3) = ss(1:3,1:3)
       end do


!      NOTE   ssinv(ii,jj)=ssinv(ii,jj)+gprimd(ii,kk)*rprimd(jj,ll)*symrec(ll,kk,isym)

       do isppol=1,nsppol
         
!        for each tensor component, rotate the cartesian directions of phonon modes
         do itensor = 1, 9
!          multiply by the ss matrices
           tmp_mat2(:,:,:) = zero
           tmp_mat(:,:,:) = reshape(gkk_qpt_tmp(:,itensor,:,isppol),&
&           (/2,nbranch,nbranch/))
           call DGEMM ('N','N',nbranch,nbranch,nbranch,&
&           one,ss_allatoms,nbranch,tmp_mat(1,:,:),nbranch,zero,&
&           tmp_mat2(1,:,:),nbranch)
           call DGEMM ('N','N',nbranch,nbranch,nbranch,&
&           one,ss_allatoms,nbranch,tmp_mat(2,:,:),nbranch,zero,&
&           tmp_mat2(2,:,:),nbranch)

           call DGEMM ('N','T',nbranch,nbranch,nbranch,&
&           one,tmp_mat2(1,:,:),nbranch,ss_allatoms,nbranch,zero,&
&           tmp_mat(1,:,:),nbranch)
           call DGEMM ('N','T',nbranch,nbranch,nbranch,&
&           one,tmp_mat2(2,:,:),nbranch,ss_allatoms,nbranch,zero,&
&           tmp_mat(2,:,:),nbranch)

           gkk_qpt_tmp(:,itensor,:,isppol) = reshape (tmp_mat, (/2,nbranch*nbranch/))
         end do ! itensor

!        for each cartesian direction/phonon mode, rotate the tensor components
         do imode = 1, nbranch*nbranch
           tmp_tensor2(:,:,:) = zero
           tmp_tensor(:,:,:) = reshape(gkk_qpt_tmp(:,:,imode,isppol),&
&           (/2,3,3/))
           call DGEMM ('N','N',3,3,3,&
&           one,sscart,3,tmp_tensor(1,:,:),3,zero,&
&           tmp_tensor2(1,:,:),3)
           call DGEMM ('N','T',3,3,3,&
&           one,tmp_tensor2(1,:,:),3,sscart,3,zero,&
&           tmp_tensor(1,:,:),3)

           call DGEMM ('N','N',3,3,3,&
&           one,sscart,3,tmp_tensor(2,:,:),3,zero,&
&           tmp_tensor2(2,:,:),3)
           call DGEMM ('N','T',3,3,3,&
&           one,tmp_tensor2(2,:,:),3,sscart,3,zero,&
&           tmp_tensor(2,:,:),3)

           gkk_qpt_tmp(:,:,imode,isppol) = reshape (tmp_tensor, (/2,9/)) ! modified by BX
         end do ! imode

!        add to gkk_qpt_new
         gkk_qpt_new(:,:,:,isppol) = gkk_qpt_new(:,:,:,isppol) + gkk_qpt_tmp(:,:,:,isppol)

       end do ! end isppol do

     end do ! end isym do
   end do ! end itim do

   ABI_CHECK(neqqpt>0,'no q-points found equivalent to iqpt ')

!  divide by number of equivalent qpts found
   gkk_qpt_new = gkk_qpt_new/neqqpt


!  copy the symmetrized version into all the equivalent qpoints, appropriately transformed
   do itim=1,2
     do isym=1,crystal%nsym
!      ieqqpt is sent onto iqpt by itim/isym
       ieqqpt = qpttoqpt(itim,isym,iqpt)

       if (symmetrized_qpt(ieqqpt) /= -1) cycle
       gkk_qpt_tmp = zero

!      use symrec matrices to get inverse transform from isym^{-1}
       do ii=1,3
         do jj=1,3
           sscart(ii,jj)=0.0_dp
           do kk=1,3
             do ll=1,3
!              Use inverse of symop matrix here to get back to ieqqpt (inv+transpose is in symrec and in gprimd)
!              sscart(ii,jj)=sscart(ii,jj)+rprimd(ii,kk)*symrec(ll,kk,isym)*gprimd(ll,jj)
               sscart(ii,jj)=sscart(ii,jj)+rprimd(ii,kk)*symrec(ll,kk,isym)*gprimd(ll,jj)
             end do
           end do
         end do
       end do
       if (ep_scalprod==1) then
         do ii=1,3
           do jj=1,3
             ss(ii,jj)=0.0_dp
             do kk=1,3
               do ll=1,3
!                Use inverse of symop matrix here to get back to ieqqpt (inv+transpose is in symrec and in gprimd)
                 ss(ii,jj)=ss(ii,jj)+rprimd(ii,kk)*symrec(ll,kk,isym)*gprimd(ll,jj)
               end do
             end do
           end do
         end do
       else
         do ii=1,3
           do jj=1,3
             ss(ii,jj) = symrel(jj,ii,isym)
           end do
         end do
       end if

       ss_allatoms(:,:) = zero
       do iatom=1,crystal%natom
         ancestor_iatom = crystal%indsym(4,isym,iatom)
         ss_allatoms((ancestor_iatom-1)*3+1:(ancestor_iatom-1)*3+3,&
&         (iatom-1)*3+1:          (iatom-1)*3+3) = ss(1:3,1:3)
       end do

!      ! Use inverse of symop matrix here to get back to ieqqpt
!      ssinv(ii,jj)=ssinv(ii,jj)+gprimd(ii,kk)*rprimd(jj,ll)*symrel(kk,ll,isym)

       do isppol=1,nsppol
         do itensor = 1, 9
!          multiply by the ss^{-1} matrices
           tmp_mat2(:,:,:) = zero
           tmp_mat(:,:,:) = reshape(gkk_qpt_new(:,itensor,:,isppol),&
&           (/2,nbranch,nbranch/))


           call DGEMM ('N','N',nbranch,nbranch,nbranch,&
&           one,ss_allatoms,nbranch,tmp_mat(1,:,:),nbranch,zero,&
&           tmp_mat2(1,:,:),nbranch)
           call DGEMM ('N','N',nbranch,nbranch,nbranch,&
&           one,ss_allatoms,nbranch,tmp_mat(2,:,:),nbranch,zero,&
&           tmp_mat2(2,:,:),nbranch)

           call DGEMM ('N','T',nbranch,nbranch,nbranch,&
&           one,tmp_mat2(1,:,:),nbranch,ss_allatoms,nbranch,zero,&
&           tmp_mat(1,:,:),nbranch)
           call DGEMM ('N','T',nbranch,nbranch,nbranch,&
&           one,tmp_mat2(2,:,:),nbranch,ss_allatoms,nbranch,zero,&
&           tmp_mat(2,:,:),nbranch)


           gkk_qpt_tmp(:,itensor,:,isppol) = reshape (tmp_mat, (/2,nbranch*nbranch/))
         end do ! itensor

!        for each cartesian direction/phonon mode, rotate the tensor components
         do imode = 1, nbranch*nbranch
           tmp_tensor2(:,:,:) = zero
           tmp_tensor(:,:,:) = reshape(gkk_qpt_tmp(:,:,imode,isppol),&
&           (/2,3,3/))
           call DGEMM ('N','N',3,3,3,&
&           one,sscart,3,tmp_tensor(1,:,:),3,zero,&
&           tmp_tensor2(1,:,:),3)
           call DGEMM ('N','T',3,3,3,&
&           one,tmp_tensor2(1,:,:),3,sscart,3,zero,&
&           tmp_tensor(1,:,:),3)

           call DGEMM ('N','N',3,3,3,&
&           one,sscart,3,tmp_tensor(2,:,:),3,zero,&
&           tmp_tensor2(2,:,:),3)
           call DGEMM ('N','T',3,3,3,&
&           one,tmp_tensor2(2,:,:),3,sscart,3,zero,&
&           tmp_tensor(2,:,:),3)

!          gkk_qpt_new(:,:,imode,isppol) = reshape (tmp_tensor, (/2,9/)) ! Modified by BX
           gkk_qpt_tmp(:,:,imode,isppol) = reshape (tmp_tensor, (/2,9/)) ! Modified by BX
         end do ! imode

         if (gkk_flag (1,1,isppol,ieqqpt) == -1) then
           gkk_flag (:,:,isppol,ieqqpt) = 0
         end if

       end do ! end isppol do


!      save symmetrized matrices for qpt ieqqpt
       gamma_qpt_tr(:,:,:,:,ieqqpt) = gkk_qpt_tmp(:,:,:,:)

       symmetrized_qpt(ieqqpt) = 1

     end do ! end isym do 
   end do ! end itim do

 end do
!end iqpt do

 ABI_DEALLOCATE(gkk_qpt_new)
 ABI_DEALLOCATE(gkk_qpt_tmp)

end subroutine complete_gamma_tr
!!***

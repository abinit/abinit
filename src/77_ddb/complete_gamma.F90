!{\src2tex{textfont=tt}}
!!****f* ABINIT/complete_gamma
!!
!! NAME
!! complete_gamma
!!
!! FUNCTION
!! Use the set of special q points calculated by the Monkhorst & Pack Technique.
!! Check if all the informations for the q points are present in the input gamma matrices.
!! Generate the gamma matrices (already summed over the FS) of the set of q points which
!! samples homogeneously the entire Brillouin zone.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2017 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! qpttoqpt = qpoint index mapping under symops
!!
!! OUTPUT
!! gamma_qpt = in/out: set of gamma matrix elements completed and symmetrized 
!!    gamma_qpt(2,nbranch**2,nsppol,nqpt_full)
!!
!! PARENTS
!!      elphon,integrate_gamma_alt,m_phgamma
!!
!! CHILDREN
!!      zgemm
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine complete_gamma(Cryst,nbranch,nsppol,nqptirred,nqpt_full,ep_scalprod,qirredtofull,qpttoqpt,gamma_qpt)

 use defs_basis
 use defs_elphon
 use m_profiling_abi
 use m_errors

 use m_crystal,    only : crystal_t

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'complete_gamma'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsppol,nbranch,nqptirred,nqpt_full,ep_scalprod
 type(crystal_t),intent(in) :: Cryst
!arrays
 integer,intent(in) :: qirredtofull(nqptirred)
 integer,intent(in) :: qpttoqpt(2,Cryst%nsym,nqpt_full)
 real(dp), intent(inout) :: gamma_qpt(2,nbranch**2,nsppol,nqpt_full)

!Local variables-------------------------------
!scalars
 integer :: ibranch,ieqqpt,ii,natom,nsym,iqpt,isppol,isym
 integer :: itim,jbranch,jj,kk,ll,neqqpt,iatom,ancestor_iatom,iqpt_fullbz
!arrays
 integer :: symmetrized_qpt(nqpt_full)
 integer :: gkk_flag(nbranch,nbranch,nsppol,nqpt_full)
 real(dp) :: ss(3,3)
 real(dp) :: tmp_mat(2,nbranch,nbranch)
 real(dp) :: tmp_mat2(2,nbranch,nbranch)
 real(dp) :: ss_allatoms(2,nbranch,nbranch)
 real(dp) :: c_one(2), c_zero(2)
 real(dp),allocatable :: gkk_qpt_new(:,:,:),gkk_qpt_tmp(:,:,:)

! *********************************************************************

 c_one = (/one,zero/)
 c_zero = (/zero,zero/)

 natom = Cryst%natom
 nsym  = Cryst%nsym

!Generation of the gkk matrices relative to the q points
!of the set which samples the entire Brillouin zone

!set up flags for gamma_qpt matrices we have
 gkk_flag = -1
 do iqpt=1,nqptirred
   iqpt_fullbz = qirredtofull(iqpt)
   gkk_flag(:,:,:,iqpt_fullbz) = 1
 end do

 symmetrized_qpt(:) = -1

 ABI_ALLOCATE(gkk_qpt_new,(2,nbranch**2,nsppol))
 ABI_ALLOCATE(gkk_qpt_tmp,(2,nbranch**2,nsppol))

 do iqpt=1,nqpt_full
!  
!  Already symmetrized?
   if (symmetrized_qpt(iqpt) == 1) cycle

   gkk_qpt_new(:,:,:) = zero

!  loop over qpoints equivalent to iqpt
   neqqpt=0
!  do not use time reversal symmetry to complete the qpoints:
!  do not know what happens to the gamma matrices
!  11/2011: MJV: time reversal is needed here if inversion is absent
!  - used in read_gkk and all reductions of q-points by symmetry.

   do itim=1,2
     do isym=1,nsym
!      ieqqpt is sent onto iqpt by itim/isym
       ieqqpt = qpttoqpt(itim,isym,iqpt)


       if (gkk_flag(1,1,1,ieqqpt) == -1) cycle
!      if we have information on this qpt
!      iqpt is equivalent to ieqqpt: get it from file or memory
       gkk_qpt_tmp(:,:,:) = gamma_qpt(:,:,:,ieqqpt)

       neqqpt=neqqpt+1

!      
!      MJV note 02/2010:
!      the correspondence of symrel and symrec in the different cases, symmetrizing there
!      and back, has been fixed in the cases with and without scalprod (ie cartesian
!      and reduced real space coordinates) with respect to a calculation with no symmetries
!      I believe everything is settled, but still do not know why the 2 versions of the ss
!      matrices here use different rel/rec, instead of just being multiplied by the rprim gprim...
!      
       if (ep_scalprod==1) then
         do ii=1,3
           do jj=1,3
             ss(ii,jj)=zero
             do kk=1,3
               do ll=1,3
                 ss(ii,jj)=ss(ii,jj)+Cryst%rprimd(ii,kk)*Cryst%symrel(kk,ll,isym)*Cryst%gprimd(ll,jj)
               end do
             end do
           end do
         end do
       else
         do ii=1,3
           do jj=1,3
             ss(ii,jj) = Cryst%symrec(ii,jj,isym)
           end do
         end do
       end if

       ss_allatoms(:,:,:) = zero
       do iatom=1,natom
         ancestor_iatom = Cryst%indsym(4,isym,iatom)
         ss_allatoms(1, (ancestor_iatom-1)*3+1:(ancestor_iatom-1)*3+3, (iatom-1)*3+1:(iatom-1)*3+3) = ss(1:3,1:3)
       end do


!      NOTE   ssinv(ii,jj)=ssinv(ii,jj)+Cryst%gprimd(ii,kk)*rprimd(jj,ll)*Cryst%symrec(ll,kk,isym)

       do isppol=1,nsppol
!        multiply by the ss matrices
         tmp_mat2(:,:,:) = zero
         tmp_mat(:,:,:) = reshape(gkk_qpt_tmp(:,:,isppol),(/2,nbranch,nbranch/))

         call ZGEMM ('N','N',nbranch,nbranch,nbranch,&
&         c_one,ss_allatoms,nbranch,tmp_mat,nbranch,c_zero,tmp_mat2,nbranch)

         call ZGEMM ('N','T',nbranch,nbranch,nbranch,&
&         c_one,tmp_mat2,nbranch,ss_allatoms,nbranch,c_zero,tmp_mat,nbranch)

!        add to gkk_qpt_new
         do ibranch =1,nbranch
           do jbranch =1,nbranch
             gkk_qpt_new(:,(jbranch-1)*nbranch+ibranch,isppol) = &
&             gkk_qpt_new(:,(jbranch-1)*nbranch+ibranch,isppol) + tmp_mat(:,jbranch,ibranch)
           end do
         end do
       end do ! isppol
!      
     end do ! isym
   end do ! itim
!  
   ABI_CHECK(neqqpt>0,'no q-points found equivalent to iqpt ')
!  Divide by number of equivalent qpts found.
   gkk_qpt_new(:,:,:) = gkk_qpt_new(:,:,:)/neqqpt

!  copy the symmetrized version into all the equivalent qpoints, appropriately transformed
   do itim=1,2
     do isym=1,nsym
!      ieqqpt is sent onto iqpt by itim/isym
       ieqqpt = qpttoqpt(itim,isym,iqpt)

       if (symmetrized_qpt(ieqqpt) /= -1) cycle
       gkk_qpt_tmp(:,:,:) = zero

!      use symrec matrices to get inverse transform from isym^{-1}
       if (ep_scalprod==1) then
         do ii=1,3
           do jj=1,3
             ss(ii,jj)=zero
             do kk=1,3
               do ll=1,3
!                Use inverse of symop matrix here to get back to ieqqpt (inv+transpose is in symrec and in gprimd)
                 ss(ii,jj)=ss(ii,jj)+Cryst%rprimd(ii,kk)*Cryst%symrec(ll,kk,isym)*Cryst%gprimd(ll,jj)
               end do
             end do
           end do
         end do
       else
         do ii=1,3
           do jj=1,3
             ss(ii,jj) = Cryst%symrel(jj,ii,isym)
           end do
         end do
       end if

       ss_allatoms(:,:,:) = zero
       do iatom=1,natom
         ancestor_iatom = Cryst%indsym(4,isym,iatom)
         ss_allatoms(1, (ancestor_iatom-1)*3+1:(ancestor_iatom-1)*3+3, (iatom-1)*3+1:(iatom-1)*3+3) = ss(1:3,1:3)
       end do

!      ! Use inverse of symop matrix here to get back to ieqqpt
!      ssinv(ii,jj)=ssinv(ii,jj)+gprimd(ii,kk)*rprimd(jj,ll)*Cryst%symrel(kk,ll,isym)

       do isppol=1,nsppol
!        multiply by the ss^{-1} matrices
         tmp_mat2(:,:,:) = zero
         tmp_mat(:,:,:) = reshape(gkk_qpt_new(:,:,isppol),(/2,nbranch,nbranch/))

         call ZGEMM ('N','N',nbranch,nbranch,nbranch,&
&         c_one,ss_allatoms,nbranch,tmp_mat,nbranch,c_zero,tmp_mat2,nbranch)

         call ZGEMM ('N','T',nbranch,nbranch,nbranch,&
&         c_one,tmp_mat2,nbranch,ss_allatoms,nbranch,c_zero,tmp_mat,nbranch)

!        FIXME: the following could just be a reshape
         do ibranch =1,nbranch
           do jbranch =1,nbranch
             gkk_qpt_tmp(:,(jbranch-1)*nbranch+ibranch,isppol) =&
&             tmp_mat(:,jbranch,ibranch)
           end do
         end do
         if (gkk_flag (1,1,isppol,ieqqpt) == -1) gkk_flag (:,:,isppol,ieqqpt) = 0
       end do ! end isppol do

!      save symmetrized matrices for qpt ieqqpt
       gamma_qpt(:,:,:,ieqqpt) = gkk_qpt_tmp(:,:,:)

       symmetrized_qpt(ieqqpt) = 1

     end do !isym
   end do !itim
 end do !iqpt

 ABI_DEALLOCATE(gkk_qpt_new)
 ABI_DEALLOCATE(gkk_qpt_tmp)

end subroutine complete_gamma
!!***

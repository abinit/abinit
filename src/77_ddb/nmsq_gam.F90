!{\src2tex{textfont=tt}}
!!****f* ABINIT/nmsq_gam
!!
!! NAME
!! nmsq_gam
!!
!! FUNCTION
!!  Calculate gamma matrices keeping full dependence on bands
!!  from original h1_mat_el_sq matrix elements (no averaging over
!!  bands near the Fermi surface)
!!
!! COPYRIGHT
!! Copyright (C) 2004-2018 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   displ_red = phonon mode displacement vectors, post-multiplied by gprim matrix
!!     (ie. turned to reduced coordinates)
!!   eigvec = phonon eigenvectors
!!   elph_ds = datastructure with gkk matrix elements
!!   FSfullpqtofull = mapping of k+q to k
!!   kpt_phon = coordinates of kpoints near to FS
!!   h1_mat_el_sq = matrix elements $<psi_{k+q,m} | H^{1} | psi_{k,n}>$ squared
!!   iqptirred = index of present qpoint
!!
!! OUTPUT
!!   accum_mat = matrix for accumulating FS average of gkk (gamma matrix -> linewidths)
!!   accum_mat2 = matrix for accumulating FS average of gamma matrix with good prefactors
!!
!! PARENTS
!!      normsq_gkq
!!
!! CHILDREN
!!      gam_mult_displ,zgemm
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine nmsq_gam (accum_mat,accum_mat2,displ_red,eigvec,elph_ds,FSfullpqtofull,&
&  h1_mat_el_sq,iqptirred)

 use defs_basis
 use defs_elphon
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nmsq_gam'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iqptirred
 type(elph_type),intent(inout) :: elph_ds
!arrays
 integer,intent(in) :: FSfullpqtofull(elph_ds%k_phon%nkpt,elph_ds%nqpt_full)
 real(dp),intent(in) :: displ_red(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp),intent(in) :: eigvec(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp),intent(inout) :: &
& h1_mat_el_sq(2,elph_ds%nFSband*elph_ds%nFSband,elph_ds%nbranch*elph_ds%nbranch,elph_ds%k_phon%my_nkpt,elph_ds%nsppol)
 real(dp),intent(inout) :: accum_mat(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)
 real(dp),intent(inout) :: accum_mat2(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)

!Local variables-------------------------------
! tmp variables for diagonalization
!scalars
 integer :: ikpt_phon,ikpt_phonq,ib1,ib2,ibeff,ibranch,isppol,ipert1
 integer :: jbranch
 integer :: iqpt_fullbz
 integer :: ik_this_proc
 real(dp) :: sd1,sd2
 character(len=500) :: message
!arrays
 real(dp) :: gkq_1band(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: tmp_mat2(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: zgemm_tmp_mat(2,elph_ds%nbranch,elph_ds%nbranch)

! *************************************************************************

 if (elph_ds%ep_keepbands == 0) then
   write (message,'(a,i0)')' elph_ds%ep_keepbands should be 1 while is ',elph_ds%ep_keepbands
   MSG_ERROR(message)
 end if

!MG20060603 NOTE:
!accum_mat and accum_mat2 are real, the imaginary part is used for debugging purpose
!accum_mat2 is used to store the phonon-linewidhts before interpolation

 iqpt_fullbz = elph_ds%qirredtofull(iqptirred)
 write(std_out,*) 'nmsq_gam : iqptirred = ', iqptirred

 do isppol=1,elph_ds%nsppol
   do ik_this_proc =1, elph_ds%k_phon%my_nkpt
     ikpt_phon = elph_ds%k_phon%my_ikpt(ik_this_proc)

     ikpt_phonq = FSfullpqtofull(ikpt_phon,iqpt_fullbz)

     do ib1=1,elph_ds%nFSband
       sd1 = elph_ds%k_phon%wtk(ib1,ikpt_phon,isppol) !weights for distance from the fermi surface

       do ib2=1,elph_ds%nFSband
         sd2 = elph_ds%k_phon%wtk(ib2,ikpt_phonq,isppol) !weights for distance from the fermi surface
         ibeff = ib2+elph_ds%nFSband*(ib1-1)

         gkq_1band(:,:,:) = zero

         zgemm_tmp_mat= reshape (h1_mat_el_sq(:,ibeff,:,ik_this_proc,isppol),(/2,elph_ds%nbranch,elph_ds%nbranch/))

         call gam_mult_displ(elph_ds%nbranch, displ_red, zgemm_tmp_mat, tmp_mat2)

!        sum over bands 
         do ipert1=1,elph_ds%nbranch
           gkq_1band(1,ipert1,ipert1) = gkq_1band(1,ipert1,ipert1) + tmp_mat2(1,ipert1,ipert1)
         end do

!        summing over k points and bands, still diagonal in jbranch
         accum_mat(:,:,:,isppol) = accum_mat(:,:,:,isppol) + gkq_1band(:,:,:)*sd1*sd2

!        MG20060603 : summing over bands and kpoints with weights to calculate the phonon linewidth
         do jbranch=1,elph_ds%nbranch
           accum_mat2(:,jbranch,jbranch,isppol) = accum_mat2(:,jbranch,jbranch,isppol) + gkq_1band(:,jbranch,jbranch)*sd1*sd2
         end do
!        END MG


!        now turn to cartesian coordinates

!        Final Gamma matrix (hermitian) = E * D_g * E^{+}
!        Where E^{+} is the hermitian conjugate of the eigenvector matrix E
!        And D_g is the diagonal matrix of values of gamma for this qpoint

!        Here gkq_1band is indexed with real phonon modes (not atom+idir)
!        turn gkq_1band to atom+cartesian coordinates (instead of normal coordinates for qpoint)
         tmp_mat2(:,:,:) = zero
         do ibranch =1,elph_ds%nbranch
           do jbranch =1,elph_ds%nbranch
             tmp_mat2(1,ibranch,jbranch) = tmp_mat2(1,ibranch,jbranch) + &
&             eigvec(1,ibranch,jbranch) * gkq_1band(1,jbranch,jbranch)
             tmp_mat2(2,ibranch,jbranch) = tmp_mat2(2,ibranch,jbranch) + &
&             eigvec(2,ibranch,jbranch) * gkq_1band(1,jbranch,jbranch)
           end do
         end do
         gkq_1band(:,:,:) = zero

!        here eigvec is transposed and complexconjugated.
         zgemm_tmp_mat=zero
         call zgemm('n','c',elph_ds%nbranch,elph_ds%nbranch,elph_ds%nbranch,cone,&
&         tmp_mat2,elph_ds%nbranch,eigvec,elph_ds%nbranch,czero,zgemm_tmp_mat,elph_ds%nbranch)

         gkq_1band = zgemm_tmp_mat

!        gamma matrix contribution in cartesian coordinates (ie interpolatable form)
         h1_mat_el_sq(:,ibeff,:,ik_this_proc,isppol) = reshape(gkq_1band,(/2,elph_ds%nbranch*elph_ds%nbranch/))

       end do
     end do
!    END loop over bands ib1 ib2

   end do
!  END loop over kpt_phon
 end do
!END loop over nsppol


end subroutine nmsq_gam
!!***

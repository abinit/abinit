!{\src2tex{textfont=tt}}
!!****f* ABINIT/nmsq_pure_gkk
!!
!! NAME
!! nmsq_pure_gkk
!!
!! FUNCTION
!!  Calculate gamma matrices for pure gkk case, ie when the
!!  scalar product with the displacement vector is done later
!!  Sum over bands is carried out later.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2018 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   displ_red = phonon displacement in reduced coordinates (used to calculate the ph linewidth)
!!   elph_ds = datastructure with gkk matrix elements
!!   FSfullpqtofull = mapping of k+q to k
!!   kpt_phon = coordinates of kpoints near to FS
!!   h1_mat_el_sq = matrix elements $<psi_{k+q,m} | H^{1} | psi_{k,n}>$ squared
!!   iqptirred = index of present qpoint
!!
!! OUTPUT
!!   elph_ds%gkq filled
!!   accum_mat = matrix for accumulating FS average of gkk (gamma matrix -> linewidths)
!!   accum_mat2 = complex array whose real part contains the phonon linewidth
!!
!! PARENTS
!!      normsq_gkq
!!
!! CHILDREN
!!      gam_mult_displ
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine nmsq_pure_gkk(accum_mat,accum_mat2,displ_red,elph_ds,FSfullpqtofull,&
&   h1_mat_el_sq,iqptirred)

 use defs_basis
 use defs_elphon
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nmsq_pure_gkk'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iqptirred
 type(elph_type),intent(inout) :: elph_ds
!arrays
 integer,intent(in) :: FSfullpqtofull(elph_ds%k_phon%nkpt,elph_ds%nqpt_full)
 real(dp),intent(in) :: displ_red(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp),intent(inout) :: &
& h1_mat_el_sq(2,elph_ds%nFSband*elph_ds%nFSband,elph_ds%nbranch*elph_ds%nbranch,elph_ds%k_phon%my_nkpt,elph_ds%nsppol)
 real(dp),intent(inout) :: accum_mat(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)
 real(dp),intent(inout) :: accum_mat2(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)

!Local variables-------------------------------
!scalars
 integer :: ikpt_phon,ikpt_phonq,ib1,ib2,ibeff,ipert1,isppol
 integer :: iqpt_fullbz
 integer :: ik_this_proc
 real(dp) :: sd1,sd2
 character(len=500) :: message
!arrays
 real(dp) :: gkq_sum_bands(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: tmp_mat2(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: zgemm_tmp_mat(2,elph_ds%nbranch,elph_ds%nbranch)

! *************************************************************************

 if (elph_ds%ep_keepbands /= 1) then
   message = ' elph_ds%ep_keepbands should be 1 to keep bands!'
   MSG_ERROR(message)
 end if

 iqpt_fullbz = elph_ds%qirredtofull(iqptirred)

!h1_mat_el_sq is already fine here - nothing to do


!MG20060603 NOTE:
!accum_mat and accum_mat2 are real, the imaginary part is used for debugging purpose
!accum_mat2 is used to store the phonon-linewidhts before interpolation

!MJV 20070525 NOTE:
!in some of the nmsq routines, in particular this one, the work done to
!calculate accum_mat,accum_mat2 is completely superfluous and will be re-done
!on the interpolated values.
!MG uses them for the QPT output, however, so keep it for consistency for the
!moment.

 do isppol=1,elph_ds%nsppol
   do ik_this_proc =1, elph_ds%k_phon%my_nkpt
     ikpt_phon = elph_ds%k_phon%my_ikpt(ik_this_proc)

     ikpt_phonq = FSfullpqtofull(ikpt_phon,iqpt_fullbz)

     gkq_sum_bands(:,:,:) = zero

!    gkq_sum_bands = \sum_{ib1,ib2} \langle k+q \mid H^{(1)}_{q,\tau_i,\alpha_i} \mid k   \rangle
!    \cdot \langle k   \mid H^{(1)}_{q,\tau_j,\alpha_j} \mid k+q \rangle
!    where ibranch -> \tau_i,\alpha_i  and  jbranch -> \tau_j,\alpha_j

     do ib1=1,elph_ds%nFSband

       sd1 = elph_ds%k_phon%wtk(ib1,ikpt_phon,isppol)      !  weights for distance from the fermi surface

       do ib2=1,elph_ds%nFSband

         sd2 = elph_ds%k_phon%wtk(ib2,ikpt_phonq,isppol)  !  weights for distance from the fermi surface
         ibeff = ib2+(ib1-1)*elph_ds%nFSband

         gkq_sum_bands = gkq_sum_bands + &
&         sd1*sd2*reshape(h1_mat_el_sq(:,ibeff,:,ik_this_proc,isppol),(/2,elph_ds%nbranch,elph_ds%nbranch/))

       end do !ib2
     end do !ib1
!    END loops over bands


     accum_mat(:,:,:,isppol) = accum_mat(:,:,:,isppol) + gkq_sum_bands(:,:,:)
   end do
!  END loop over kpt_phon

!  MG20060603
!  do scalar product with the displ_red to calculate the ph lwdth before interpolation (stored in accum_mat2)

   zgemm_tmp_mat = accum_mat(:,:,:,isppol)

   call gam_mult_displ(elph_ds%nbranch, displ_red, zgemm_tmp_mat, tmp_mat2)

   do ipert1=1,elph_ds%nbranch
     accum_mat2(1,ipert1,ipert1,isppol) = accum_mat2(1,ipert1,ipert1,isppol) + tmp_mat2(1,ipert1,ipert1)
   end do

!  ENDMG

 end do ! isppol

end subroutine nmsq_pure_gkk
!!***

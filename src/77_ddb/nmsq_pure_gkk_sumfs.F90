!{\src2tex{textfont=tt}}
!!****f* ABINIT/nmsq_pure_gkk_sumfs
!!
!! NAME
!! nmsq_pure_gkk_sumfs
!!
!! FUNCTION
!!  Calculate gamma matrices for pure gkk case, i.e, when the
!!  scalar product with the displacement vector is done later
!!  Sum over bands is carried out now.
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

subroutine nmsq_pure_gkk_sumfs(accum_mat,accum_mat2,displ_red,elph_ds,FSfullpqtofull,h1_mat_el_sq,iqptirred)

 use defs_basis
 use defs_elphon
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nmsq_pure_gkk_sumfs'
 use interfaces_77_ddb, except_this_one => nmsq_pure_gkk_sumfs
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iqptirred
 type(elph_type),intent(in) :: elph_ds
!arrays
 integer,intent(in) :: FSfullpqtofull(elph_ds%k_phon%nkpt,elph_ds%nqpt_full)
 real(dp),intent(in) :: displ_red(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp),intent(inout) :: &
& h1_mat_el_sq(2,elph_ds%nFSband*elph_ds%nFSband,elph_ds%nbranch*elph_ds%nbranch,elph_ds%k_phon%my_nkpt,elph_ds%nsppol)
 real(dp),intent(inout) :: accum_mat(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)
 real(dp),intent(inout) :: accum_mat2(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)

!Local variables-------------------------------
!scalars
 integer :: ikpt_phon,ikpt_phonq,ib1,ib2,ibeff,ipert1,isppol,iqpt_fullbz
 integer :: nbranch,nsppol,nFSband,nkpt_phon
 integer :: ik_this_proc
 real(dp) :: sd1,sd2
 !character(len=500) :: message
!arrays
 real(dp) :: gkq_sum_bands(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: tmp_mat2(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: zgemm_tmp_mat(2,elph_ds%nbranch,elph_ds%nbranch)

! *************************************************************************

 if (elph_ds%ep_keepbands /= 0) then
   MSG_BUG('ep_keepbands should be 0 to average over bands!')
 end if

 nbranch   = elph_ds%nbranch
 nsppol    = elph_ds%nsppol
 nFSband   = elph_ds%nFSband
 nkpt_phon = elph_ds%k_phon%nkpt

 iqpt_fullbz = elph_ds%qirredtofull(iqptirred)

!MG20060603 NOTE:
!accum_mat and accum_mat2 are real, the imaginary part is used for debugging purpose
!accum_mat2 is used to store the phonon-linewidhts before interpolation

 do isppol=1,nsppol
   do ik_this_proc =1, elph_ds%k_phon%my_nkpt
     ikpt_phon = elph_ds%k_phon%my_ikpt(ik_this_proc)

!    
!    The index of k+q in the BZ.
     ikpt_phonq = FSfullpqtofull(ikpt_phon,iqpt_fullbz)
!    
!    gkq_sum_bands = 
!    \sum_{ib1,ib2} <k+q| H^{(1)}_{q,\tau_i,\alpha_i} |k> \cdot <k| H^{(1)}_{q,\tau_j,\alpha_j}|k+q>
!    
!    where ibranch = (\tau_i,\alpha_i) and  jbranch = (\tau_j,\alpha_j).
     gkq_sum_bands(:,:,:) = zero

     do ib1=1,nFSband
       sd1 = elph_ds%k_phon%wtk(ib1,ikpt_phon,isppol)      !  weights for distance from the fermi surface

       do ib2=1,nFSband
         sd2 = elph_ds%k_phon%wtk(ib2,ikpt_phonq,isppol)  !  weights for distance from the fermi surface
         ibeff=ib2+(ib1-1)*nFSband

         gkq_sum_bands = gkq_sum_bands + &
&         sd1*sd2* reshape(h1_mat_el_sq(:,ibeff,:,ik_this_proc,isppol),(/2,nbranch,nbranch/))
       end do !ib2
     end do !ib1
!    
!    gamma matrix contribution in reduced coordinates (ie interpolatable form)
!    The sum over Fermi surface bands is done here, and fed into (ib1,ib2)=(1,1)
     h1_mat_el_sq(:,1,:,ik_this_proc,isppol) = reshape(gkq_sum_bands,(/2,nbranch**2/))

     accum_mat(:,:,:,isppol) = accum_mat(:,:,:,isppol) + gkq_sum_bands(:,:,:)
   end do ! kpt_phon
 end do ! isppol
!
!MG20060603
!do scalar product wit displ_red to calculate the ph lwdth before interpolation (stored in accum_mat2)
 do isppol=1,nsppol
   zgemm_tmp_mat = accum_mat(:,:,:,isppol)
!  
   call gam_mult_displ(nbranch, displ_red, zgemm_tmp_mat, tmp_mat2)

   do ipert1=1,nbranch
     accum_mat2(1,ipert1,ipert1,isppol) = accum_mat2(1,ipert1,ipert1,isppol) + tmp_mat2(1,ipert1,ipert1)
   end do
!  
 end do

end subroutine nmsq_pure_gkk_sumfs
!!***

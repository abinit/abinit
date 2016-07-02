!!****f* ABINIT/gam_mult_displ
!!
!! NAME
!! gam_mult_displ
!!
!! FUNCTION
!! This routine takes the bare gamma matrices and multiplies them
!!  by the displ_red matrices (related to the scalprod variable)
!!
!! COPYRIGHT
!! Copyright (C) 2009-2016 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   nbranch = number of phonon branches (3*natom)
!!   displ_red = phonon mode displacement vectors in reduced coordinates.
!!   gam_bare = bare gamma matrices before multiplication
!!
!! OUTPUT
!!   gam_now = output gamma matrices multiplied by displacement matrices
!!
!! PARENTS
!!      get_tau_k,m_phgamma,mka2f,mka2f_tr,mka2f_tr_lova,mkph_linwid,nmsq_gam
!!      nmsq_gam_sumfs,nmsq_pure_gkk,nmsq_pure_gkk_sumfs,normsq_gkq
!!
!! CHILDREN
!!      zgemm
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine gam_mult_displ(nbranch, displ_red, gam_bare, gam_now)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gam_mult_displ'
!End of the abilint section

 implicit none

!Arguments -------------------------------
 integer, intent(in)  :: nbranch
 real(dp), intent(in)  :: displ_red(2,nbranch,nbranch)
 real(dp), intent(in)  :: gam_bare(2,nbranch,nbranch)
 real(dp), intent(out) :: gam_now(2,nbranch,nbranch)

!Local variables -------------------------
 real(dp) :: zgemm_tmp_mat(2,nbranch,nbranch)

! *********************************************************************

 gam_now = zero

 call zgemm('c','n',nbranch,nbranch,nbranch,cone,&
& displ_red,nbranch,gam_bare,&
& nbranch,czero,zgemm_tmp_mat,nbranch)

 call zgemm('n','n',nbranch,nbranch,nbranch,cone,&
& zgemm_tmp_mat,nbranch,displ_red,&
& nbranch,czero,gam_now,nbranch)

end subroutine gam_mult_displ
!!***

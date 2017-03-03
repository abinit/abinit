!!****f* ABINIT/k_neighbors
!!
!! NAME
!!   k_neighbors
!!
!! FUNCTION
!!   find 8 neighbors of given k-point on a coarse grid, and return
!!   them along with relative k-shift within coarse grid cell
!!
!! COPYRIGHT
!! Copyright (C) 2010-2017 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   kpt        = k-point to be interpolated to, in full BZ
!!   kptrlatt   = lattice vectors for coarse k-grid
!!   invrankkpt = rank list to find k-points
!!
!! OUTPUT
!!   rel_kpt = k-point coordinates renormalized to coarse grid cell
!!   kpt_phon_indices = indices of k-points on corners of cell
!!
!! PARENTS
!!      integrate_gamma_alt
!!
!! CHILDREN
!!      get_rank_1kpt,interpol3d_indices,wrap2_zero_one
!!
!! SOURCE 

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine k_neighbors (kpt, kptrlatt,kptrank_t, rel_kpt, kpt_phon_indices)

 use defs_basis
 use m_kptrank
 use m_profiling_abi

 use m_numeric_tools,  only : wrap2_zero_one, interpol3d_indices

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'k_neighbors'
!End of the abilint section

 implicit none

! inputs
 real(dp), intent(in) :: kpt(3)
 integer, intent(in) :: kptrlatt(3,3)
 type(kptrank_type), intent(in) :: kptrank_t

! outputs
 real(dp), intent(out) :: rel_kpt(3)
 integer, intent(out) :: kpt_phon_indices(8)

! local vars
 integer :: symrankkpt
 integer :: ir1,ir2,ir3, pr1,pr2,pr3
 real(dp) :: redkpt(3), cornerkpt(3), res

! *************************************************************************

!wrap fine kpt to [0,1]
 call wrap2_zero_one(kpt(1),redkpt(1),res)
 call wrap2_zero_one(kpt(2),redkpt(2),res)
 call wrap2_zero_one(kpt(3),redkpt(3),res)
!find 8 indices of points neighboring ikpt_phon, for interpolation
 call interpol3d_indices (redkpt,kptrlatt(1,1),kptrlatt(2,2),kptrlatt(3,3), &
& ir1,ir2,ir3, pr1,pr2,pr3)

!transpose ir pr to ikpt_phon indices
!order of kpt_phons:
!ir1 ir2 ir3
 cornerkpt = (/real(ir1-1)/kptrlatt(1,1),real(ir2-1)/kptrlatt(2,2), real(ir3-1)/kptrlatt(3,3)/)
 call get_rank_1kpt (cornerkpt,symrankkpt,kptrank_t)
 kpt_phon_indices(1) = kptrank_t%invrank(symrankkpt)
!pr1 ir2 ir3 
 cornerkpt = (/real(pr1-1)/kptrlatt(1,1),real(ir2-1)/kptrlatt(2,2), real(ir3-1)/kptrlatt(3,3)/)
 call get_rank_1kpt (cornerkpt,symrankkpt,kptrank_t)
 kpt_phon_indices(2) = kptrank_t%invrank(symrankkpt)
!ir1 pr2 ir3 
 cornerkpt = (/real(ir1-1)/kptrlatt(1,1),real(pr2-1)/kptrlatt(2,2), real(ir3-1)/kptrlatt(3,3)/)
 call get_rank_1kpt (cornerkpt,symrankkpt,kptrank_t)
 kpt_phon_indices(3) = kptrank_t%invrank(symrankkpt)
!pr1 pr2 ir3 
 cornerkpt = (/real(pr1-1)/kptrlatt(1,1),real(pr2-1)/kptrlatt(2,2), real(ir3-1)/kptrlatt(3,3)/)
 call get_rank_1kpt (cornerkpt,symrankkpt,kptrank_t)
 kpt_phon_indices(4) = kptrank_t%invrank(symrankkpt)
!ir1 ir2 pr3 
 cornerkpt = (/real(ir1-1)/kptrlatt(1,1),real(ir2-1)/kptrlatt(2,2), real(pr3-1)/kptrlatt(3,3)/)
 call get_rank_1kpt (cornerkpt,symrankkpt,kptrank_t)
 kpt_phon_indices(5) = kptrank_t%invrank(symrankkpt)
!pr1 ir2 pr3 
 cornerkpt = (/real(pr1-1)/kptrlatt(1,1),real(ir2-1)/kptrlatt(2,2), real(pr3-1)/kptrlatt(3,3)/)
 call get_rank_1kpt (cornerkpt,symrankkpt,kptrank_t)
 kpt_phon_indices(6) = kptrank_t%invrank(symrankkpt)
!ir1 pr2 pr3 
 cornerkpt = (/real(ir1-1)/kptrlatt(1,1),real(pr2-1)/kptrlatt(2,2), real(pr3-1)/kptrlatt(3,3)/)
 call get_rank_1kpt (cornerkpt,symrankkpt,kptrank_t)
 kpt_phon_indices(7) = kptrank_t%invrank(symrankkpt)
!pr1 pr2 pr3 
 cornerkpt = (/real(pr1-1)/kptrlatt(1,1),real(pr2-1)/kptrlatt(2,2), real(pr3-1)/kptrlatt(3,3)/)
 call get_rank_1kpt (cornerkpt,symrankkpt,kptrank_t)
 kpt_phon_indices(8) = kptrank_t%invrank(symrankkpt)

!retrieve the gkq matrix for all q, at the neighbor k vectors
 rel_kpt(1) = redkpt(1)*kptrlatt(1,1)-real(ir1-1)
 rel_kpt(2) = redkpt(2)*kptrlatt(2,2)-real(ir2-1)
 rel_kpt(3) = redkpt(3)*kptrlatt(3,3)-real(ir3-1)

end subroutine k_neighbors
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/expibi
!! NAME
!! expibi
!!
!! FUNCTION
!! Routine that computes exp(i (-b_ket).R) at each site.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2018 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  dkvecs(3) :: $\Delta k$ increment
!!  natom :: number of atoms in unit cell
!!  xred(natom,3) :: reduced coordinates of atoms in unit cell
!!
!! OUTPUT
!!  calc_expibi(2,natom) :: phase factors at each atom for vector shift
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      initberry,overlap_k1k2_paw
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine expibi(calc_expibi,dkvecs,natom,xred)

 use m_profiling_abi
 use m_errors
 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'expibi'
!End of the abilint section

 implicit none

!Arguments---------------------------
!scalars
 integer,intent(in) :: natom
 real(dp),intent(out) :: calc_expibi(2,natom)
!arrays
 real(dp),intent(in) :: dkvecs(3),xred(3,natom)

!Local variables---------------------------
!scalars
 integer :: iatom
 real(dp) :: bdotr

! *************************************************************************

 calc_expibi(:,:) = zero

 !calc_expibi(2,natom)
 !used for PAW field calculations (distributed over atomic sites)
 !stores the on-site phase factors arising from
 !$\langle\phi_{i,k}|\phi_{j,k+\sigma_k k_k}\rangle$
 !where $\sigma = \pm 1$. These overlaps arise in various Berry
 !phase calculations of electric and magnetic polarization. The on-site
 !phase factor is $\exp[-i\sigma_k k_k)\cdot I]$ where
 !$I$ is the nuclear position. 

 do iatom = 1, natom

    !    note the definition used for the k-dependence of the PAW basis functions:
    !$|\phi_{i,k}\rangle = exp(-i k\cdot r)|\phi_i\rangle
    !    see Umari, Gonze, and Pasquarello, PRB 69,235102 Eq. 23. 
    bdotr = DOT_PRODUCT(xred(1:3,iatom),-dkvecs(1:3))
    !    here is exp(i b.R) for the given site
    calc_expibi(1,iatom) = cos(two_pi*bdotr)
    calc_expibi(2,iatom) = sin(two_pi*bdotr)

 end do ! end loop on natom

end subroutine expibi
!!***

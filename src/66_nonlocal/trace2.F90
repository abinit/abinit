!{\src2tex{textfont=tt}}
!!****f* ABINIT/trace2
!! NAME
!! trace2
!!
!! FUNCTION
!! Sum indices to compute trace of rank 2 tensor gxa related to l=2
!! $trace=sum_{i,j} {gxa(i,j) gmet(i,j)}$
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  gxa(2,6) = $sum_{G} {e^(2 \pi i G cdot t} {{f_2(|k+G|)} \over {|k+G|^2}} (k+G) cdot (k+G) C(G_{nk})}$
!!  gmet(3,3)=(symmetric) metric tensor in reciprocal space (bohr^-2)
!!
!! OUTPUT
!!  trace(2)=sum_{i,j} {gxa(i,j) gmet(i,j)}$ (Re and Im).
!!
!! NOTES
!! Here index 6 refers to vector components
!! of (k+G) but note tensor is symmetric=>only 6 components.
!! The components are given in the order 11 22 33 32 31 21.
!! The initial 2 handles the Re and Im parts.
!!
!! PARENTS
!!      nonlop_pl
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine trace2(gxa,gmet,trace)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'trace2'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: gmet(3,3),gxa(2,6)
 real(dp),intent(out) :: trace(2)

!Local variables-------------------------------
!scalars
 integer :: reim

! *************************************************************************

!Write out index summation, Re and Im parts
 do reim=1,2
   trace(reim)=gxa(reim,1)*gmet(1,1)+gxa(reim,2)*gmet(2,2)+&
&   gxa(reim,3)*gmet(3,3)+&
&   2.0d0*(gxa(reim,4)*gmet(3,2)+gxa(reim,5)*gmet(3,1)+&
&   gxa(reim,6)*gmet(2,1))
 end do

end subroutine trace2
!!***
